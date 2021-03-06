#' Access the Tonsil Atlas data (RNA, ATAC, CITE, Spatial)
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
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom HDF5Array HDF5Array
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment assay
#'
#' @examples
HCATonsilData <- function(assayType = "RNA", cellType = "All", processedCounts = TRUE) {
  # Sanity checks
  allowedAssays <- c("RNA", "ATAC", "CITE", "Spatial")
  if (!(assayType %in% allowedAssays)) {
    stop(
      "'assayType' must be included in ",
      paste(allowedAssays, collapse = ", "))
  }
  allowedTypes <- listCellTypes(assayType = assayType)
  if (!(cellType %in% allowedTypes)) {
    stop(
      "'cellType' must be included in ",
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

  # Update annotations
  sce <- updateAnnotation(
    sce = sce,
    refAnnotation = "20220215",
    newAnnotation = "20220619"
  )

  # We did some last minute changes to NBC-MBC before publication, let's
  # reannotate and change the UMAP coords to map with the manuscript
  if (cellType == "NBC-MBC") {
    sce <- sce[, NBC_MBC_annotation_df$barcode]
    sce$annotation_20220619 <- NBC_MBC_annotation_df$annotation_20220619
    umap_df <- as.matrix(NBC_MBC_annotation_df[, c("UMAP_1", "UMAP_2")])
    reducedDim(sce, "UMAP") <- umap_df
  }

  # Similarly, let us update NBC/MBC annotation in case cellType = All
  if (cellType == "All") {
    annot <- sce$annotation_20220619
    names(annot) <- colnames(sce)
    annot[NBC_MBC_annotation_df$barcode] <- NBC_MBC_annotation_df$annotation_20220619
    sce$annotation_20220619 <- annot
  }

  # Return
  sce
}

#' Information on the Tonsil Data
#'
#' Get information on the individual datasets of the Tonsil Atlas data
#'
#' @inheritParams HCATonsilData
#'
#' @return An info message is generated, and its content is returned invisibly.
#'
#' @author Federico Marini

#' @export
#'
#' @examples
#' HCATonsilDataInfo(assayType = "RNA", cellType = "epithelial")
HCATonsilDataInfo <- function(assayType = "RNA", cellType = "All") {
  # Sanity checks
  allowedAssays <- c("RNA", "ATAC", "CITE", "Spatial")
  if (!(assayType %in% allowedAssays)) {
    stop(
      "'assayType' must be included in ",
      paste(allowedAssays, collapse = ", "))
  }
  allowedTypes <- listCellTypes(assayType = assayType)
  if (!(cellType %in% allowedTypes)) {
    stop(
      "'cellType' must be included in ",
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

  # Fetch info instead of creating full sce object
  rd_sce <- eh[eh$rdatapath == filePaths["rowData"]][[1]]
  cd_sce <- eh[eh$rdatapath == filePaths["colData"]][[1]]

  # buildup the information string
  info_sce <- sprintf(
    "Dataset: %s (assay selected: %s)\nMeasurements available for %d genes in %d cells",
    cellType,
    assayType,
    nrow(rd_sce),
    nrow(cd_sce)
  )

  message(info_sce)

  return(invisible(info_sce))
}
