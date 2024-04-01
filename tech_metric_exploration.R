######## Example ######## 

# TODO: Make a wrapper function for all QC

# Set up metadata table for samples
# Two samples for testing
#sample_meta <- data.frame(
#  sample_id = c("Xen_PCa_9224", "CosMx_PCa_9224"),
#  path = c(
#    "/Volumes/Charmander/spatial_data/new_data/Xenium_Prostate_9224_Rep1/",
#    "/Volumes/Charmander/spatial_data/new_data/CosMx_Prostate_9224_Rep1/"
#  ),
#  platform = c("Xenium", "CosMx")
#)

# Full directory of samples
path <- "/Volumes/Charmander/touchstone/"
sample_meta <- data.frame(
  sample_id = list.files(path),
  path = paste0(path, list.files(path)),
  platform = stringr::word(list.files(path), 3,3, sep="_")
)
sample_meta$platform[1:4] <- rep("CR", 4)
sample_meta$platform[sample_meta$platform == "XR"] <- "Xenium"
sample_meta$platform[sample_meta$platform == "CR"] <- "CosMx"

# Load objects from table
obj_list <- list()
for(i in 1:nrow(sample_meta)){
  obj_list[[i]] <- readSpatial(sample_id = sample_meta[i, "sample_id"],
                               path = sample_meta[i, "path"],
                               platform = sample_meta[i, "platform"])
}


## Annotate
#ref loaded from side script of snPATHO-seq from matched sample--will clean
#obj_list <- lapply(obj_list, annotateData, ref = ref)

## QC
# Probe set size
panel_size <- data.frame(
  sample_id = sample_meta$sample_id,
  platform = sample_meta$platform,
  value = unlist(lapply(obj_list, nrow))
)

# Total cells
cell_count <- data.frame(
  sample_id = sample_meta$sample_id,
  platform = sample_meta$platform,
  value = unlist(lapply(obj_list, ncol))
)

# Transcripts per cell
tx_per_cell <- do.call(rbind, lapply(obj_list, getTxPerCell))

# Transcripts per um2
tx_per_um2 <- do.call(rbind, lapply(obj_list, getTxPerArea))

# Transcripts per nucleus
tx_per_nuc <- do.call(rbind, lapply(obj_list, getTxPerNuc))

# Transcript per cell per 100 genes in panel
tx_per_cell_norm <- tx_per_cell
tx_per_cell_norm$value <- round(tx_per_cell_norm$value / panel_size$value,
                                digits = 2)

# Transcripts per cell with intersecting genes
common_genes <- Reduce(intersect, lapply(obj_list, rownames))
tx_per_cell_intersect <- do.call(rbind, lapply(obj_list, getTxPerCell, features=common_genes))

# Fraction transcripts in cells
tx_fraction_in_cell <- do.call(rbind, lapply(obj_list, getCellTxFraction))

# Mean signal-noise ratio
signal_ratio <- do.call(rbind, lapply(obj_list, getMeanSignalRatio))

# Expression
mean_expression <- do.call(rbind, lapply(obj_list, getMeanExpression))

# Max expression distribution
max_exp <- do.call(rbind, lapply(obj_list, getMaxDetection))

# MECR
mecr <- do.call(rbind, lapply(obj_list, getMECR))

# Moran's I
morans <- do.call(rbind, lapply(obj_list, getMorans))

# Silhouette Width
#silhouette <- do.call(rbind, lapply(obj_list, getSilhouetteWidth))

# PLOT
p0 <- plotSampleLabel(sample_meta)
p1 <- plotPanelSize(panel_size)
p2 <- plotCellCount(cell_count)
p3 <- plotTxPerCell(tx_per_cell)
p4 <- plotTxPerArea(tx_per_um2) 
p5 <- plotTxPerNuc(tx_per_nuc)
p6 <- plotTxPerCellNorm(tx_per_cell_norm) 
p65 <- plotTxPerCell(tx_per_cell_intersect) + xlab("Per cell\ncommon\ngenes")
p7 <- plotFractionTxInCell(tx_fraction_in_cell)
p8 <- plotSignalRatio(signal_ratio) 
p9 <- plotMeanExpression(mean_expression)
p10 <- plotMECR(mecr)
p11 <- plotMorans(morans)

p <- cowplot::plot_grid(p0, p1, p2, p3, p4, p5, p6, p65, p7, p8, p9, p10, p11, 
                        ncol=13, align='h',
                        rel_widths = c(1, 0.4, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.75, 0.75, 0.8, 0.3, 0.8, 0.8))
cowplot::save_plot(p, filename="~/Downloads/benchmark_qc_small.pdf",
                   base_width=20, base_height=20)

