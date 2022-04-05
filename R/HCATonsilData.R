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
  allowedAssays <- c("RNA", "ATAC", "CITE", "Spatial")
  if (!(assayType %in% allowedAssays)) {
    stop(
      "'assay_type' must be included in ",
      paste(allowedAssays, collapse = ", "))
  }
  allowedTypes <- listCellTypes(assayType = assayType)
  if (!(cellType %in% allowedTypes)) {
    stop(
      "'assay_type' must be included in ",
      paste(allowedTypes, collapse = ", "))
  }
  eh <- ExperimentHub::ExperimentHub()
  host <- file.path("HCATonsilData", "1.0", assayType)
}
