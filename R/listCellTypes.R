#' List available cell types for a dataset of the tonsil cell atlas
#'
#' @param assayType Either 'RNA', 'ATAC', 'CITE', or 'Spatial'
#'
#' @export
#'
#' @return A character vector with the available cell types for the indicated
#'   dataset.
#' @examples
#' listCellTypes(assayType = "RNA")
#'
listCellTypes <- function(assayType) {
  if (assayType == "RNA") {
    c(
      "All", "NBC-MBC", "GCBC", "PC","CD4-T", "Th", "CD8-T", "ILC-NK",
      "myeloid", "FDC", "epithelial", "PDC", "preB", "preT"
    )
  } else if (assayType == "ATAC") {
    c("All", "NBC/MBC", "GCBC", "PC","CD4 T", "CD8 T", "ILC/NK")
  } else if (assayType == "CITE") {
    c("All", "CD4 T")
  } else if (assayType == "Spatial") {
    "All"
  } else {
    stop("Invalid 'dataset' (must be either 'RNA', 'ATAC', 'CITE' or 'Spatial'")
  }
}
