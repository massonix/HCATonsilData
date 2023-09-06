#' Update cell type annotation
#'
#' Annotations are dynamic by nature. As more eyes look into the data and
#' newer literature comes out, we expect annotations to be refined over time.
#' We have accommodated this by allowing us and users to add newer annotations
#' to the SingleCellExperiment objects.
#'
#' @param sce A SingleCellExperiment object obtained using HCATonsilData
#'   function.
#' @param refAnnotation string specifying the date of the annotation to use as
#'   reference.
#' @param newAnnotation string specifying the suffix to add to the new column
#'   (annotation_*).
#'
#' @return A SingleCellExperiment object with additional columns
#'   (annotation_*) that contains more annotations.
#' @export
#'
#' @examples
#' # TODO
updateAnnotation <- function(sce,
                             refAnnotation = "20220215",
                             newAnnotation = "20220619") {
  dict_20220619 <- c(
    "GC-Tfh-0X40" = "GC-Tfh-OX40",
    "non-GC-Tf-regs" = "Eff-Tregs-IL32",
    "GC-Tf-regs" = "Tfr",
    "IFN CD8 T" = "IFN+ CD8 T",
    "CXCR6+ RM CD8 T" = "RM CD8 activated T",
    "IL7R MMP12 macrophages" = "MMP Slancytes",
    "C1Q HLA macrophages" = "C1Q Slancytes",
    "SELENOP FUCA1 PTGDS macrophages" = "SELENOP Slancytes",
    "ITGAX ZEB2 macrophages" = "ITGAX Slancytes",
    "Mast cells" = "Mast",
    "LZ FDC" = "FDC",
    "DZ FDC" = "COL27A1+ FDC",
    "DZ_Sphase" = "DZ early Sphase",
    "DZ_Sphase_HistoneHigh" = "DZ late Sphase",
    "DZ_G2M_HistoneHigh" = "DZ early G2Mphase",
    "DZ_G2M_CCNBHigh"= "DZ late G2Mphase",
    "DZ-cell cycle exit" = "DZ cell cycle exit",
    "DZ-nonproliferative" = "DZ cell cycle exit",
    "DZ-nonproliferative_FOXP1hi"= "DZ non proliferative",
    "DZ/LZ" = "DZ_LZ transition",
    "DZ/LZ" = "DZ_LZ transition",
    "LZ" = "LZ",
    "LZ-BCL2A1 neg"= "LZ",
    "LZ-DZ-re-entry early commitment" = "LZ_DZ reentry commitment",
    "LZ-proliferative_BCL2A1pos" = "LZ proliferative",
    "LZ-proliferative_BCL2A1neg" = "LZ_DZ transition",
    "MBC-like_nonproli" = "Precursor MBCs",
    "MBC-like_proli1" = "Precursor MBCs",
    "MBC-like_proli2"= "Reactivated proliferative MBCs",
    "MBC-like_proli3" = "Reactivated proliferative MBCs",
    "MBC-like_FCRL4+"= "Reactivated proliferative MBCs",
    "PC-precursors" = "PC committed Light Zone GCBC",
    "class switch MBC" = "csMBC",
    "Neutrophil Granulocytes" = "Neutrophils"
  )
  oldAnnot <- paste("annotation", refAnnotation, sep = "_")
  newAnnot <- paste("annotation", newAnnotation, sep = "_")
  sce[[newAnnot]] <- as.character(sce[[oldAnnot]])
  dictSub <- dict_20220619[names(dict_20220619) %in% unique(sce[[oldAnnot]])]
  for (cellType in names(dictSub)) {
    sce[[newAnnot]][sce[[newAnnot]] == cellType] <- dictSub[cellType]
  }
  sce
}


#' TODO
#'
#' TODO
#'
#' @details TODO
#'
#' @name NBC_MBC_annotation_df
#' @docType data
NULL
