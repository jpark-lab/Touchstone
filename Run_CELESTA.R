library(CELESTA)
library(Rmixmod)
library(spdep)
library(ggplot2)
library(reshape2)
library(zeallot)
library(dplyr)
library(Seurat)

# using CELESTA to annotate protein data
# from https://github.com/plevritis-lab/CELESTA

annotate_PRT_insitu_data <- function(Sobj, ref) {
  ## load ref profile matrix from project yaml
  prior_marker_info <- read.csv(project$paths$ref_celesta_directory)
  query_mat <- Sobj[["RNA"]]$counts # marked RNA because protein counts are also stored as "RNA" assay
  imaging_data <- t(query_mat) %>% as.data.frame()
  imaging_data$X <- Sobj$x_slide_mm
  imaging_data$Y <- Sobj$y_slide_mm
  overlapping_markers <- colnames(prior_marker_info)[colnames(prior_marker_info) %in% colnames(imaging_data)]
  CelestaObj <- CreateCelestaObject(project_title = "CosMx_protein", prior_marker_info[,c("celltype", "Lineage_level", overlapping_markers)], imaging_data[,c("X", "Y", overlapping_markers)])
  CelestaObj <- FilterCells(CelestaObj, high_marker_threshold = 0.9, low_marker_threshold = 0.4)
  CelestaObj <- AssignCells(CelestaObj, max_iteration = 10, cell_change_threshold = 0.01,
                          high_expression_threshold_anchor = high_marker_threshold_anchor,
                          low_expression_threshold_anchor = low_marker_threshold_anchor,
                          high_expression_threshold_index = high_marker_threshold_iteration,
                          low_expression_threshold_index = low_marker_threshold_iteration)
  return(CelestaObj)
}