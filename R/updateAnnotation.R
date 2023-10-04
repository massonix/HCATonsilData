#' Update cell type annotation
#'
#' Annotations are dynamic by nature. As more experts look into the data and
#' newer literature comes out, we expect annotations to be refined over time.
#' We have accommodated this by allowing us and users to add new annotations
#' to the SingleCellExperiment objects. If you want to propose a new
#' annotation based on your experience or new evidence, please open an issue
#' at https://github.com/massonix/HCATonsilData/issues.
#'
#' @param sce A SingleCellExperiment object obtained using HCATonsilData
#'   function.
#' @param refAnnotation string specifying the date of the annotation to use as
#'   reference.
#' @param newAnnotation string specifying the suffix to add to the new column
#'   (annotation_*).
#'
#' @return A SingleCellExperiment object with additional an additional column
#'   (annotation_*) that contains more annotations.
#' @export
#'
#' @examples
#' update the annotation from preprint (version 1.0) to publication (version 2.0)
#' sce <- updateAnnotation(
#'   sce,
#'   refAnnotation = "20220619",
#'    newAnnotation = "20230508"
#' )

updateAnnotation <- function(sce,
                             refAnnotation = "20220215",
                             newAnnotation = "20220619") {
  oldAnnot <- paste("annotation", refAnnotation, sep = "_")
  newAnnot <- paste("annotation", newAnnotation, sep = "_")
  selectedDict <- paste0("dict_", refAnnotation, "_to_", newAnnotation)
  selectedDict <- HCATonsilData::annotations_dictionary[[selectedDict]]
  sce[[newAnnot]] <- as.character(sce[[oldAnnot]])
  selectedDict <- selectedDict[names(selectedDict) %in% unique(sce[[oldAnnot]])]
  for (cellType in names(selectedDict)) {
    sce[[newAnnot]][sce[[newAnnot]] == cellType] <- selectedDict[cellType]
  }
  sce
}
