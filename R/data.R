#' Annotation dictionary for naive and memory B cells (NBC-MBC)
#'
#' In version 1 of the atlas, we changed the clusters and annotation for NBC-MBC
#' prior to publishing the preprint. Since, there is not a 1:1 mapping between
#' clusters, here we provide the correspondence between the annotations in
#' February 2022 and the preprint (June 2022) for each NBC-MBC cell barcode.
#'
#' This data is used by the updateAnnotation function in B cells.
#'
#'
#' @format A dataframe with 112478 NBC-MBC and 5 variables
#' \describe{
#'   \item{barcode}{cell barcode}
#'   \item{names_level_5_clusters_eta}{Cell annotation given by the B-cell annotation team in February 2022}
#'   \item{annotation_20220619}{Cell annotation given by the B-cell annotation team June 2022 (preprint)}
#'   \item{UMAP_1}{UMAP1 coordinates}
#'   \item{UMAP_2}{UMAP2 coordinates}
#' }
#'
#' @source \url{https://zenodo.org/record/8373756}
#' @name NBC_MBC_annotation_df
#' @docType data
globalVariables("NBC_MBC_annotation_df")



#' Annotation dictionary to keep track of annotations changes and versions
#' 
#' Annotations are opinionated and dynamic by nature. They are subjected to 
#' change as more experts look into the data and more evidence is published.
#' HCATonsilData has been designed to account for that. We timestamp annotations
#' and plan to include additional annotations if users propose them and we
#' validate them.
#'
#'
#' @format A list of named vectors with the correspondence between time-stamped
#'   annotations
#'
#' @source This file was created with the script in
#' "inst/scripts/make-annotation-dictionary.R
#' @name annotation_dictionary
#' @docType data
globalVariables("annotation_dictionary")


#' Donor metadata
#' 
#' Data frame that contains all metadata information for the 17 donors included
#' in the tonsil atlas.
#' 
#'
#' @format A dataframe with 17 observations and 8 variables
#'
#' @source Check the tonsil atlas publication
#' @name donor_metadata
#' @docType data
globalVariables("donor_metadata")


#' Color palettes
#' 
#' Color palettes used in the paper for all cell types
#' 
#' 
#' @format A list with the HEX color codes for each cell type
#' 
#' @source This file was created with the script in
#' "inst/scripts/make-annotation-colors.R
#' @name colors_20230508
#' @docType data
globalVariables("colors_20230508")