#' Title Access the Tonsil Atlas data (RNA, ATAC, CITE, Spatial)
#'
#' The data was downloaded from Zenodo
#' https://zenodo.org/record/6340174#.YkwQwX9BxH4
#'
#' @param assayType One of "RNA", "ATAC", "CITE" or "Spatial".
#' @param cellType A character vector with the cell type available.
#'    A list of available cell types can be obtained using
#'    \code{listCellTypes(assay_type)}.
#' @param processedCounts Logical scalar. If \code{TRUE}, include the processed
#'   counts in addition to the raw counts in the SingleCellExperiment object.
#'
#' @return A named list of
#'   \linkS4class{SingleCellExperiment} objects (one per cell_type requested
#'   via \code{cell_type})
#'
#' @author Ramon Massoni-Badosa

#' @export
#'
#' @examples
HCATonsilData <- function(assayType, cellType, processedCounts = TRUE) {
  # Sanity checks
  allowedAssays <- c("RNA", "ATAC", "CITE", "Spatial")
  if (!(assayType %in% allowedAssays)) {
    stop(
      "'assay_type' must be included in ",
      paste(allowedAssays, collapse = ", "))
  }
  allowedTypes <- listCellTypes(assayType = assayType)
  if (!(cellType %in% allowedTypes)) {
    stop(
      "'assayType' must be included in ",
      paste(allowedTypes, collapse = ", "))
  }

  # Initialize ExperimentHub object
  eh <- ExperimentHub::ExperimentHub()
  host <- file.path("HCATonsilData", "1.0", assayType)

  # Define paths
  suffixes <- c("counts.h5", "processed.h5", "rowdata.rds", "coldata.rds",
                "pca.rds", "harmony.rds", "umap.rds")
  names(suffixes) <- c("counts", "processedCounts", "rowData", "colData",
                       "pca", "harmony", "umap")
  filePaths <- sapply(suffixes, \(.) {
    x <- file.path(host, paste(cellType, assayType, ., sep = "_"))
    x
  })
  for (x in filePaths) {
    if (sum(x == eh$rdatapath) > 1) {
      stop("Input matched more than one entry!")
    }
  }

  # Create SingleCellExperiment
  cts <- HDF5Array::HDF5Array(
    eh[eh$rdatapath == filePaths["counts"]][[1]],
    name = "counts"
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = cts),
    rowData = eh[eh$rdatapath == filePaths["rowData"]][[1]],
    colData = eh[eh$rdatapath == filePaths["colData"]][[1]]
  )
  SingleCellExperiment::reducedDims(sce) <- list(
    PCA = eh[eh$rdatapath == filePaths["pca"]][[1]],
    UMAP = eh[eh$rdatapath == filePaths["umap"]][[1]]
  )
  if (sum(filePaths["harmony"] == eh$rdatapath) == 1) {
    harmony <- eh[eh$rdatapath == filePaths["harmony"]][[1]]
    SingleCellExperiment::reducedDims(sce)[["HARMONY"]] <- harmony
  }
  if (processedCounts) {
    prccts <- HDF5Array::HDF5Array(
      eh[eh$rdatapath == filePaths["processedCounts"]][[1]],
      name = "processed"
  )
  SummarizedExperiment::assay(sce, "logcounts", withDimnames = FALSE) <- prccts
  }
  sce
}
