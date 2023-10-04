#' List available cell types for a dataset of the tonsil cell atlas
#'
#' @param assayType Either 'RNA', 'ATAC', 'Multiome', 'CITE', or 'Spatial'
#' @param version Version of the tonsil atlas data: '1.0' (preprint) or '2.0'
#'   or '2.0' (publication, default). Note that, for version 2.0, CD8-T and
#'   ILC-NK are combined in a single 'Cytotoxic' object.
#'
#' @export
#'
#' @return A character vector with the available cell types for the indicated
#'   dataset.
#' @examples
#' listCellTypes(assayType = "RNA", version = "2.0")
#'

listCellTypes <- function(assayType=c("RNA", "ATAC", "Multiome", "CITE", "Spatial"),
                          version="2.0") {
  assayType <- match.arg(assayType)
  version <- match.arg(version, c("1.0", "2.0"))
  if (assayType %in% c("ATAC", "Multiome", "CITE", "Spatial")) {
    cellTypes <- "All"
  } else if(assayType == "RNA") {
    cellTypes <- switch(
      version,
      "1.0" = c("All", "NBC-MBC", "GCBC", "PC","CD4-T", "Th", "CD8-T", "ILC-NK",
                "myeloid", "FDC", "epithelial", "PDC", "preB", "preT"),
      "2.0" = c("All", "NBC-MBC", "GCBC", "PC", "CD4-T", "Cytotoxic", "myeloid",
                "FDC", "epithelial", "PDC")
    )
  }
  cellTypes
}