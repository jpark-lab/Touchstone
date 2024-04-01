library(Seurat)
library(tidyverse)
library(ggplot2)
library(shadowtext)
library(scales)
library(cowplot)
library(data.table)
library(Matrix)
library(matrixStats)
library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(bluster)
library(NMI)
library(spacexr)
library(InSituType)
library(BiocParallel)

######## I/O ######## 

# Functions assume data is stored in specific formats:
# Xenium: Typical Xenium bundle; CosMx: Flat CSV (TAP style); MERSCOPE: Not implemented yet

# readSpatial() reads in data from either Xenium, CosMx, of MERSCOPE.
# It outputs a seurat object with some common metadata e for downstream comparison.
# Regardless of platform, data is stored in an assay named "RNA" for convenient
# function
readSpatial <- function(sample_id, path, platform){
  print(paste0("Reading: ", sample_id))
  
  if(platform == "Xenium"){
    print("Loading Xenium data")
    seu_obj <- LoadXenium(path, assay = "RNA")
    seu_obj@meta.data$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data
    
    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, "cells.csv.gz")
    )
    #Set up a few defined metadata columns
    seu_obj@meta.data$cell_area <- cell_meta$cell_area
    seu_obj@meta.data$nucleus_area <- cell_meta$nucleus_area
    seu_obj@meta.data$transcript_counts <- cell_meta$transcript_counts
    seu_obj@meta.data$negprobe_counts <- cell_meta$control_probe_counts
    
    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- GetTissueCoordinates(seu_obj)
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")
    
  } else if(platform == "CosMx"){
    print("Loading CosMx data")
    seu_obj <- LoadNanostring(path, fov="fov")
    seu_obj$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data
    
    #Fix assays to separate targeting and non-targeting probes
    ## Add negative control probe assay
    sys_probes <- grep("SystemControl", rownames(seu_obj), value=T)
    neg_probes <- grep("Negative", rownames(seu_obj), value=T)
    seu_obj[["ControlProbe"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[neg_probes,]
    )
    ## Make "Nanostring" assay
    tx_probes <- rownames(seu_obj)[!rownames(seu_obj) %in% c(sys_probes, neg_probes)]
    seu_obj[["RNA"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[tx_probes,]
    )
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj[["Nanostring"]] <- NULL
    
    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, list.files(path, pattern="*metadata_file.csv.gz"))
    )
    #It's excessive, but we'll add all metadata into the object
    seu_obj@meta.data <- cbind(seu_obj@meta.data, cell_meta)
    seu_obj@meta.data$fov <- factor(paste0("FOV", seu_obj@meta.data$fov))
    seu_obj@meta.data$cell_area <- seu_obj$Area.um2
    seu_obj@meta.data$transcript_counts <- seu_obj$nCount_RNA
    seu_obj@meta.data$negprobe_counts <- seu_obj$nCount_ControlProbe
    
    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- data.frame(
      Tissue_1 = cell_meta$CenterX_global_px,
      Tissue_2 = cell_meta$CenterY_global_px
    )
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")

    
  } else if(paltform == "Merscope"){
    print("Working on support!")
    stop()
    
  } else{
    print("Not a supported platform")
    stop()
    
  }
  
  return(seu_obj)
}

# readTxMeta() simply reads in the transcript localization/metadata table
# for each platform. This table will be used by subsequent functions
readTxMeta <- function(path, platform){
  if(platform == "Xenium"){
    df <- data.table::fread(file.path(path, "transcripts.csv.gz"))
  } else if(platform == "CosMx"){
    df <- data.table::fread(file.path(path, 
                                   list.files(path, pattern = "*tx_file.csv.gz")))
    df <- unique(df)
  } else if(platform == "Merscope"){
    print("Working on support!")
    stop()
  } else{
    print("Platform not supported")
    stop()
  }
}

######## General Utilities ######## 
getPseudobulk <- function(seu_obj, celltype_meta="cell_type") {
  
  
  celltype <- factor(seu_obj@meta.data[,celltype_meta])
  names(celltype) <- colnames(seu_obj)
  mat <- seu_obj[["RNA"]]$counts
  
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(s) {
    cells <- names(celltype)[celltype==s]
    pseudobulk <- rowMeans(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}

# Automated annotation of spatial data with single-cell references using InSituType
# Assumes `ref` is a seurat object with a "cell_type" metadata column--can be overridden
# NOTE: Currently, annotateData() will remove cells with <10 transcripts
# so be cautious of the data composition changing slightly
annotateData <- function(seu_obj, ref, celltype_meta="cell_type"){
  print("Getting pseudobulk for reference")
  ref_mat <- getPseudobulk(ref)
  
  cells_keep <- colnames(seu_obj)[colSums(seu_obj[["RNA"]]$counts) > 10]
  seu_obj <- subset(seu_obj, cells = cells_keep)
  
  query_mat <- seu_obj[["RNA"]]$counts
  
  print("Annotated spatial data")
  insitutype_res <- insitutypeML(x = t(query_mat),
                                 neg = colMeans(seu_obj[["ControlProbe"]]$counts),
                                 reference_profiles = ref_mat) 
  
  seu_obj$celltype_pred <- insitutype_res$clust
  return(seu_obj)
}


######## QC Metrics ######## 

### Transcripts per cell

getTxPerCell <- function(seu_obj, #features can be explicitly defined. Defaults to all targets
                         features=NULL){ 
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  mean_tx <- mean(colSums(seu_obj[["RNA"]]$counts[features,]))
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean_tx
  )
  return(res)
} 

### Transcripts per um2
getTxPerArea <- function(seu_obj, 
                         features=NULL){ 
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  tx_count <- colSums(seu_obj[["RNA"]]$counts[features,])
  mean_tx_norm <- mean(tx_count / seu_obj$cell_area)
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean_tx_norm
  )
  return(res)
}

### Transcripts per nucleus
getTxPerNuc <- function(seu_obj, 
                        features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)
  
  # Read Tx localization data
  tx_df <- readTxMeta(path, platform)
  
  if(platform == "Xenium"){
    tx_df <- filter(tx_df, cell_id %in% colnames(seu_obj) & 
                      overlaps_nucleus == 1 &
                      feature_name %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())
    
    
  } else if(platform == "CosMx"){
    tx_df$cell_id <- paste(tx_df$cell_ID, tx_df$fov, sep="_")
    tx_df <- tx_df %>%
      filter(cell_id %in% colnames(seu_obj) & 
               CellComp == "Nuclear" &
               target %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())
    
  } else if(platform == "Merscope"){
    print("Working on support")
    
  } else{
    print("Platform not supported")
  }
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean(tx_df$nuc_counts)
  )
  
  return(res)
}

### Per Probe Mean Expression
getMeanExpression <- function(seu_obj, 
                              features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  target_df <- data.frame(
    target = features,
    value = rowMeans(seu_obj[["RNA"]]$counts[features,]),
    type = "Gene"
  )
  
  control_df <- data.frame(
    target = rownames(seu_obj[["ControlProbe"]]$counts),
    value = rowMeans(seu_obj[["ControlProbe"]]$counts),
    type = "Control"
  )
  
  res <- rbind(target_df, control_df)
  res$platform <- unique(seu_obj$platform)
  res$sample_id <- unique(seu_obj$sample_id)
  
  return(res)
} 

### log-ratio of mean gene counts to mean neg probe counts
getMeanSignalRatio <- function(seu_obj, 
                               features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
  neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)
  
  ratio <- log10(tx_means) - log10(mean(neg_probe_means))
  ratio <- mean(ratio)
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=ratio
  )
  
  return(res)
}

### Fraction of transcripts in cells
getCellTxFraction <- function(seu_obj, 
                              features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)
  
  tx_df <- readTxMeta(path, platform)
  
  if(platform == "Xenium"){
    tx_df <- filter(tx_df, feature_name %in% features)
    total_tx_count <- nrow(tx_df)
    unassigned_tx_count <- sum(tx_df$cell_id == "UNASSIGNED")
    
    cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count
    
  } else if(platform == "CosMx"){
    tx_df <- filter(tx_df, target %in% features)
    total_tx_count <- nrow(tx_df)
    unassigned_tx_count <- sum(tx_df$CellComp == "None")
    
    cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count
    
  } else if(platform == "Merscope"){
    print("Working on support")
    
  } else{
    print("Platform not supported")
  }
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=cell_tx_fraction
  )
  
  return(res)
}

##### Dynamic Range 
# Log-ratio of highest mean exp vs. mean noise
getMaxRatio <- function(seu_obj, 
                        features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
  neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)
  
  ratio <- log10(max(tx_means)) - log10(mean(neg_probe_means))
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=ratio
  )
  
  return(res)
}

# Distribution of maximal values
getMaxDetection <- function(seu_obj, 
                            features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  max_vals <- matrixStats::rowMaxs(as.matrix(seu_obj[["RNA"]]$counts[features,]))
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=max_vals,
    gene = features
  )
  
  return(res)
}

##### Mutually Exclusive Co-expression Rate (MECR) Implementation
getMECR <- function(seu_obj) {
  #This function comes from Hartman & Satija, bioRxiv, 2024
  #We are using a custom marker table. The original publication bases it on
  #scRNA-seq from matched tissue.
  marker_df <- data.frame(
    gene = c("EPCAM", "KRT19", "KRT8", 
             "CD3E", "CD3D", "CD8A", "NKG7",
             "MS4A1", "CD79A",
             "PECAM1", "CLDN5", "VWF",
             "C1QA", "C1QB", "CD14", "FCGR3A", "ITGAX", "ITGAM",
             "PDGFRA", "DPT", "COL1A1",
             "MYH11", "ACTG2"),
    cell_type = c("Epithelial", "Epithelial", "Epithelial",
                  "T", "T", "T", "T",
                  "B", "B",
                  "Endo", "Endo", "Endo",
                  "Macro", "Macro", "Macro", "Macro", "Macro", "Macro",
                  "Fibro", "Fibro", "Fibro",
                  "Muscle", "Muscle")
  )
  rownames(marker_df) <- marker_df$gene
  
  coexp.rates <- c()
  genes <- intersect(rownames(seu_obj), rownames(marker_df))
  print(paste0("Marker count: ", length(genes)))
  if (length(genes) > 25) { genes <- sample(genes, 25) }
  mtx <- as.matrix(seu_obj[['RNA']]$counts[genes, ])
  for (g1 in genes) {
    for (g2 in genes) {
      if ((g1 != g2) && (g1 > g2) && (marker_df[g1, "cell_type"] != marker_df[g2, "cell_type"])) {
        c1 <- mtx[g1, ]
        c2 <- mtx[g2, ]
        coexp.rates <- c(
          coexp.rates,
          sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0)) # >0 too liberal of an expression threshold?
      }
    }
  }
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(mean(coexp.rates), digits=3)
  )
  
  return(res)
}

##### Distribution spatial autocorrelation
getMorans <- function(seu_obj, 
                      features=NULL){
  #Requires SingleCellExperiment, SpatialFeatureExperiment, Voyager, scater
  
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  #First run for gene-targeting probes
  print("Getting Moran's I for gene-targeting probes")
  sce <- SingleCellExperiment(list(counts=seu_obj[["RNA"]]$counts[features,]),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)
  
  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                                 dist_type = "idw", k = 20, 
                                                 style = "W")
  
  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))
  
  spatial_cor <- as.data.frame(rowData(sfe))
  
  targeting <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Gene"
  )
  
  #Now run for control probes
  print("Getting Moran's I for non-targeting probes")
  sce <- SingleCellExperiment(list(counts=seu_obj[["ControlProbe"]]$counts),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)
  
  #Nearest neighbor
  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                                 dist_type = "idw", k = 20, 
                                                 style = "W")
  #Moran's I
  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))
  
  spatial_cor <- as.data.frame(rowData(sfe))
  
  control <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Control"
  )
  
  res <- rbind(targeting, control)
  
  return(res)
}

##### Cluster evaluation: silhouette width
getSilhouetteWidth <- function(seu_obj){
  print("Clustering data")
  seu_obj <- seu_obj %>%
    NormalizeData() %>%
    ScaleData()
  VariableFeatures(seu_obj) <- rownames(seu_obj)
  seu_obj <- seu_obj %>%
    RunPCA(verbose=F) %>%
    FindNeighbors(dims=1:10) %>%
    FindClusters(resolution=0.2)
  
  silhouette <- bluster::approxSilhouette(
    Embeddings(seu_obj, 'pca')[,1:10],
    clusters = seu_obj$seurat_clusters
  )
  
  silhouette <- as.data.frame(silhouette)
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(mean(silhouette$width), digits=3)
  )
  
}

getRCTD <- function(seu_obj, ref){
  # Prep reference for RCTD
  counts <- ref[["RNA"]]$counts
  cluster <- ref$cell_type
  cluster <- factor(cluster)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)
  
  # Prep obj
  # NOTE: If tx counts are low, it will fail. 
  # I've found that subseting to >10 counts tends to work. 
  # Should decide if that's done here or in a separate function
  print("Filtering any cells with <10 tx counts")
  seu_obj <- subset(seu_obj, nCount_RNA > 10)
  counts <- seu_obj[["RNA"]]$counts
  coords <- GetTissueCoordinates(seu_obj)
  coords <- filter(coords, cell %in% colnames(seu_obj)) %>%
    column_to_rownames(var = "cell")
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  # Run RCTD in full mode
  print("Running RCTD")
  RCTD <- create.RCTD(query, reference, max_cores = 6,
                      UMI_min = 0, UMI_max = Inf, counts_MIN = 0,
                      UMI_min_sigma = 50)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  RCTD@results$weights <- RCTD@results$weights / rowSums(RCTD@results$weights)
  return(RCTD)
  #seu_obj$max_weight <- rowMaxs(RCTD@results$weights)
  
}

getMaxRCTD <- function(seu_obj, ref){
  # Prep reference for RCTD
  counts <- ref[["RNA"]]$counts
  cluster <- ref$cell_type
  cluster <- factor(cluster)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)
  
  # Prep obj
  # NOTE: If tx counts are low, it will fail. 
  # I've found that subseting to >10 counts tends to work. 
  # Should decide if that's done here or in a separate function
  print("Filtering any cells with <10 tx counts")
  seu_obj <- subset(seu_obj, nCount_RNA > 10)
  counts <- seu_obj[["RNA"]]$counts
  coords <- GetTissueCoordinates(seu_obj)
  coords <- filter(coords, cell %in% colnames(seu_obj)) %>%
    column_to_rownames(var = "cell")
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  # Run RCTD in full mode
  print("Running RCTD")
  RCTD <- create.RCTD(query, reference, max_cores = 6,
                      UMI_min = 0, UMI_max = Inf, counts_MIN = 0,
                      UMI_min_sigma = 50)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  RCTD@results$weights <- RCTD@results$weights / rowSums(RCTD@results$weights)
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(rowMaxs(RCTD@results$weights), digits=3)
  )
  
  return(res)
}

#Returns data frame with expression levels for both objects for plotting
getCorrelationExp <- function(seu_obj, ref){
  common_genes <- intersect(rownames(seu_obj), rownames(ref))
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    gene = common_genes,
    value_sample = rowMeans(seu_obj[["RNA"]]$counts[common_genes,]),
    value_ref = rowMeans(ref[["RNA"]]$counts[common_genes,]),
    neg_probe_signal = mean(seu_obj[["ControlProbe"]]$counts)
  )
  return(res)
}

#Returns spearman correlation for two objects
getCorrelation <- function(seu_obj, ref){
  common_genes <- intersect(rownames(seu_obj), rownames(ref))
  cor_res <- cor(
    rowMeans(seu_obj[["RNA"]]$counts[common_genes,]),
    rowMeans(ref[["RNA"]]$counts[common_genes,]),
    method='spearman'
  )
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=cor_res
  )
  return(res)
}

# Gets correlation per cell type
# Assumes seu_obj has `celltype_pred` and ref has `cell_type`
getCellTypeCor <- function(seu_obj, ref){
  cell_types <- levels(ref$cell_type)
  
  cor_list <- list()
  for(i in 1:length(cell_types)){
    print(paste0("Correlating: ", cell_types[i]))
    cor_list[[i]] <- getCorrelation(
      subset(seu_obj, celltype_pred == cell_types[i]),
      subset(ref, cell_type == cell_types[i])
    )
  }
  cor_list <- do.call(rbind, cor_list)
  cor_list$cell_type <- cell_types
  return(cor_list)
}

# Get proportion of cells of a given cell type
# Assumes seu_obj has `celltype_pred` and ref has `cell_type`
getCellTypeProportion <- function(seu_obj, ref){
  sample_prop <- seu_obj@meta.data %>%
    mutate(cell_type = seu_obj@meta.data[,"celltype_pred"]) %>%
    group_by(cell_type) %>%
    summarize(count = n()) %>%
    mutate(prop = count / sum(count))
  ref_prop <- ref@meta.data %>%
    mutate(cell_type = ref@meta.data[,"cell_type"]) %>%
    group_by(cell_type) %>%
    summarize(count = n()) %>%
    mutate(prop = count / sum(count))
  sample_prop$group <- "In Situ"
  ref_prop$group <- "snPATHO"
  df <- bind_rows(sample_prop, ref_prop)
  df$sample_id <- unique(seu_obj$sample_id)
  df$platform <- unique(seu_obj$platform)
  return(df)
}

getClusterMetrics <- function(seu_obj, metadata_col){
  print("Clustering data")
  seu_obj <- seu_obj %>%
    NormalizeData() %>%
    ScaleData() 
  
  VariableFeatures(seu_obj) <- rownames(seu_obj)
  
  seu_obj <- seu_obj %>%
    RunPCA(verbose=F) %>%
    FindNeighbors(dims=1:20) %>%
    FindClusters(resolution=0.5) 
  
  ari <- bluster::pairwiseRand(seu_obj$seurat_clusters,
                               seu_obj$celltype_pred,
                               mode=('index'),
                               adjusted=TRUE)
  
  nmi <- NMI(
    data.frame(cellid = colnames(seu_obj), cluster = seu_obj$seurat_clusters),
    data.frame(cellid = colnames(seu_obj), cluster = seu_obj$celltype_pred)
  )
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=c(ari, nmi$value),
    metric = c("ARI", "NMI")
  )
  
}

######## Plotting ######## 
# All plots assume input is a tidy data frame with the following columns:
# 1) sample_id
# 2) platform
# 3) value (based on what is being plotted--from functions above)

plotSampleLabel <- function(sample_meta){
  df <- data.frame(
    samples = sample_meta$sample_id,
    platform = sample_meta$platform
  )
  
  p <- ggplot(df, aes(x="", y=samples)) +
    geom_text(aes(label=samples, color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) + 
    theme_void() + theme(legend.position='none')
  return(p)
}

plotPanelSize <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("Panel size") + ylab("") +
    scale_size(range = c(6,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
  
}

plotCellCount <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_text(aes(label=scales::comma(value), color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) + 
    scale_x_discrete(position='top') +
    xlab("Cell count") + ylab("") +
    theme_classic() + 
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  return(p)
}

plotTxPerCell <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("Per cell") + ylab("") +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerArea <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\num^2") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerNuc <- function(df){
  df$column <- ""
  p <- ggplot(df, aes(x=column, y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per\nnucleus") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerCellNorm <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\nper target") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotFractionTxInCell <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90', stroke=1) +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Fraction Tx in Cells") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0),
                       breaks=c(0, 0.5, 1),
                       labels = c(0, 0.5, 1),
                       limits=c(0,1)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerCell_Intersect <- function()

plotSignalRatio <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Mean log10-ratio\nexpression over noise") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotMeanExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black", 
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n detection per cell") + ylab("") +
    scale_x_log10(position='top', expand = c(0,0),
                  labels = label_log(digits = 2)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotMaxExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=platform)) +
    geom_boxplot(color="black", fill="lightgrey",
                 alpha=0.5, outlier.size=0, outlier.colour = NA) +
    scale_colour_manual(values = c("#59C134", "#14B3E6")) +
    xlab("Maximal gene\ndetection per cell") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0),
                  limits = c(0,50), oob=squish) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plotMECR <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9, "YlOrRd"),
                         limits=c(0,0.1)) +
    xlab("MECR") + ylab("") +
    scale_size(range = c(7,12), limits = c(0, 0.1)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotMorans <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black", 
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n spatial autocorrelation\n(Moran's I)") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotSilhouette <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -0.001) +
    xlab("Mean\nsilhouette width\n(Louvain res=0.2)") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotCorrelation <- function(df){
  #facet_wrap order is opposite to axis text orders, so we'll flip the levels
  df$sample_id <- factor(df$sample_id)
  df$sample_id <- factor(df$sample_id, levels = rev(levels(df$sample_id)))
  p <- ggplot(df, aes(x=value_ref, y=value_sample)) +
    geom_point(size=0.5, shape=16, stroke=0, alpha=0.75) +
    geom_abline(intercept = 0, slope = 1, linetype=2) +
    scale_x_log10(labels = label_log(digits = 2), limits=c(1e-4, 10), 
                  oob=squish, position='top') +
    scale_y_log10(labels = label_log(digits = 2), limits=c(1e-4, 10), oob=squish) +
    xlab("snPATHO-seq\ncorrelation") + ylab("") +
    facet_wrap(~sample_id, ncol=1) +
    theme_bw() +
    theme(axis.text = element_text(size=10, color="black"),
          legend.position="none",
          axis.title.x = element_text(size=12),
          strip.text = element_blank(),
          strip.background = element_blank())
  return(p)
}

plotCellTypeCor <- function(df){
  p <- ggplot(df, aes(x=value, y = sample_id)) +
    geom_jitter(shape=21, color='black', aes(fill=cell_type),
                alpha=0.5, size=3, height=0.1) +
    stat_summary(fun = mean, geom = "crossbar", width=0.5) +
    scale_x_continuous(position='top', limits=c(0, 1),
                       expand=c(0,0), breaks=c(0, 0.5, 1)) +
    xlab("Cell type\ncorrelation") + ylab("") +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 6.5, 5.5, 5.5, unit='pt')
    )
  return(p)
}

plotCellTypeProportion <- function(df){
  p <- ggplot(df, aes(x=prop, y=group)) +
    geom_bar(position="stack", stat="identity",
             color='black', alpha=0.75,
             aes(fill=cell_type), width=0.7) +
    scale_x_continuous(position='top',
                       limits=c(0,1), expand = c(0,0),
                       breaks=c(0, 0.5, 1)) +
    ylab("") + xlab("Cell type\nproportion") +
    facet_wrap(~sample_id, ncol=1) +
    theme_classic() +
    theme(
      panel.spacing = grid::unit(0, "lines"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position='none',
      axis.text.y = element_text(size=10, color="black"),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 6.5, 5.5, 5.5,
                           unit = "pt")
    )
  return(p)
}


plotARI <- function(cluster_metrics){
  df <- cluster_metrics %>% filter(metric == "ARI")
  
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=round(value, digits=2))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("ARI") + ylab("") +
    scale_size(range = c(6,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotNMI <- function(cluster_metrics){
  df <- cluster_metrics %>% filter(metric == "NMI")
  
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=round(value, digits=2))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("NMI") + ylab("") +
    scale_size(range = c(6,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotRCTD <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_boxplot(color="black", 
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 fill="lightgrey", width=0.5) +
    xlab("Max decomposition\nweight") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0),
                       limits=c(0, 1), breaks=c(0, 0.5, 1),
                       oob=squish) +
    #scale_x_log10(position='top', expand = c(0,0),
    #              labels = label_log(digits = 2)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 6.5, 5.5, 5.5, unit='pt')
    )
  return(p)
}


