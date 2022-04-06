# .makeSingleCellExperiment <- function(host, cell_type) {
#   cts <- rhdf5::h5read(path_to_cts, name = "counts") # This will need to be queried from ExperimentHub
#   mat <- rhdf5::h5read(path_to_processed, name = "processed") # This will need to be queried from ExperimentHub
#   row_data <- readRDS(path_to_rowdata)
#   col_data <- readRDS(path_to_coldata)
#   pca <- readRDS(path_to_pca)
#   harmony <- readRDS(path_to_harmony)
#   umap <- readRDS(path_to_umap)
#   sce <- SingleCellExperiment::SingleCellExperiment(
#     assays = list(counts = cts),
#     colData = coldata,
#     rowData = rowdata
#   )
#   SingleCellExperiment::reducedDims(sce) <- list(
#     PCA = pca,
#     HARMONY = harmony,
#     UMAP = umap
#   )
#   sce
# }
#
#
