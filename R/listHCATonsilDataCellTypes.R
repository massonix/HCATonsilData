#' List available cell types for a dataset of the tonsil cell atlas
#'
#' @param dataset Either 'RNA', 'ATAC', 'CITE', or 'Spatial'
#'
#' @export
#'
#' @return A character vector with the available cell types for the indicated
#'   dataset.
#' @examples
#' listHCATonsilDataCellTypes(dataset = "RNA")
#'
listHCATonsilDataCellTypes <- function(dataset) {
  if (dataset == "RNA") {
    c(
      "All", "NBC/MBC", "GCBC", "PC","CD4 T", "Th", "CD8 T", "ILC/NK",
      "myeloid", "FDC", "epithelial", "PDC", "preB", "preT"
    )
  } else if (dataset == "ATAC") {
    c("All", "NBC/MBC", "GCBC", "PC","CD4 T", "CD8 T", "ILC/NK")
  } else if (dataset == "CITE") {
    c("All", "CD4 T")
  } else if (dataset == "Spatial") {
    "All"
  } else {
    stop("Invalid 'dataset' (must be either 'RNA', 'ATAC', 'CITE' or 'Spatial'")
  }
}
