#' TonsilData_glossary
#'
#' Convenience function to read directly in the file provided as `extdata`
#'
#' @return Data frame containing the info on the cell types included in the
#' TonsilDataAtlas
#' @export
#'
#' @examples
#' glossary_df <- TonsilData_glossary()
TonsilData_glossary <- function() {
  glossary_location <- system.file("extdata", "sloglo_tabular.csv", package = "HCATonsilData")

  glossary_df <- read.table(glossary_location)

  return(glossary_df)
}




#' TonsilData_cellinfo
#'
#' @param cell_type String character, used to define the cell type for which the
#' information will be displayed. Defaults to `NULL` - if left unspecified, the
#' function returns a list of the possible options
#'
#' @return Invisible `NULL` - the information is displayed as a message in the
#' console.
#'
#' @export
#'
#' @examples
#' TonsilData_cellinfo()
#'
#' TonsilData_cellinfo("PDC")
TonsilData_cellinfo <- function(cell_type = NULL) {

  glossary_df <- TonsilData_glossary()
  slo_celltypes <- rownames(glossary_df)

  if (is.null(cell_type)) {
    message("Please select one of the following: ",
            paste0(slo_celltypes, collapse = "|"))
  } else {

    if (cell_type %in% slo_celltypes) {

      cell_msg <- paste0(
        glossary_df[cell_type, "annotation_detailed"],
        "\n------------------------------",
        "\nAnnotation Level 1: ",
        glossary_df[cell_type, "annotation_level_1"],
        "\nCell Markers: ",
        glossary_df[cell_type, "description"],
        "\nCell Markers: ",
        glossary_df[cell_type, "marker_genes"],
        "\nRelated references: ",
        glossary_df[cell_type, "related_refs"] # ,
        # "\nCell Ontology terms: ",
        # glossary_df[matched_celltype, "related_cellontology"]
      )

      message(cell_msg)


    } else {
      message("Cell type not found! Please select one of the following: ",
              paste0(slo_celltypes, collapse = "|"))
    }
  }

  invisible(NULL)
}


