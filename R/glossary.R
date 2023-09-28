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
#' head(glossary_df)
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
#' @importFrom utils browseURL read.table
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



.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"

library("htmltools")

.link_marker <- function(gene_id) {
  sprintf(
    '<a href = "http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
    gene_id,
    .actionbutton_biocstyle,
    gene_id
  )
}

.link_cellontology <- function(cl_id) {
  # https://www.ebi.ac.uk/ols/ontologies/cl/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FCL_0000084
  # actually this one is better: http://purl.obolibrary.org/obo/CL_0000000 (and works)

  cl_term <- gsub("CL:", "CL_", cl_id)

  sprintf(
    '<a href = "http://purl.obolibrary.org/obo/%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
    cl_term,
    .actionbutton_biocstyle,
    cl_id
  )
}

.link_refs <- function(refname, ref_doi) {
  sprintf(
    '<a href = "https://doi.org/%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
    ref_doi,
    .actionbutton_biocstyle,
    refname
  )
}








#' TonsilData_cellinfo_html
#'
#' @param cell_type String character, used to define the cell type for which the
#' information will be displayed. Defaults to `NULL` - if left unspecified, the
#' function returns a list of the possible options
#' @param display_plot Logical value, defines whether to include or not a plot
#' for the UMAP of the selected cell type
#' @param output_to Character value, one of "single_page" or "html_to_embed".
#' Defines in which form the information is returned, either as an individual page
#' or simply as HTML code to directly embed into other documents.
#'
#' @return Either the HTML code generated to be embedded, or a standalone HTML
#' page is created - and the location to this temp file is returned as a character
#' (default).
#' @export
#'
#' @importFrom htmltools tags
#' @importFrom utils browseURL read.table
#' @importFrom rmarkdown render
#' @importFrom base64enc dataURI
#'
#' @examples
#' TonsilData_cellinfo_html("PDC")
#' TonsilData_cellinfo_html("preB")
#' TonsilData_cellinfo_html("preT")
#' TonsilData_cellinfo_html("preT", output_to = "html_to_embed")
TonsilData_cellinfo_html <- function(cell_type = NULL,
                                     display_plot = TRUE,
                                     output_to = c("single_page", "html_to_embed")) {


  output_to <- match.arg(output_to, c("single_page", "html_to_embed"))

  glossary_df <- TonsilData_glossary()
  slo_celltypes <- rownames(glossary_df)

  if (is.null(cell_type)) {
    message("Please select one of the following: ",
            paste0(slo_celltypes, collapse = "|"))
  } else {

    if (cell_type %in% slo_celltypes) {

      cell_html_celltype <- paste0(
        tags$b(cell_type),
        tags$br(),
        tags$b("Annotation Level 1: "),
        glossary_df[cell_type, "annotation_level_1"],
        tags$hr(),
        tags$br()
      )

      cell_html_celldescription <- paste0(
        tags$b("Cell type description: "),
        glossary_df[cell_type, "description"],
        tags$br(), tags$br()
      )

      markers <- as.character(glossary_df[cell_type, "marker_genes"])
      markers_split <- unlist(strsplit(markers, split = ",", fixed = TRUE))
      markers_content <- paste(
        unlist(lapply(markers_split, function(id) .link_marker(id))),
        collapse = " ")

      cell_html_cellmarkers <- paste0(
        tags$b("Marker genes for this cell type: "),
        tags$br(),
        markers_content,
        tags$br(),
        tags$br()
      )


      refs <- as.character(glossary_df[cell_type, "related_refs"])
      refs_split <- unlist(strsplit(refs, split = ";", fixed = TRUE))
      refs_split_names <- unlist(lapply(
        strsplit(refs_split, split = "|", fixed = TRUE),
        function(arg) arg[1]))
      refs_split_dois <- unlist(lapply(
        strsplit(refs_split, split = "|", fixed = TRUE),
        function(arg) arg[2]))

      refs_content <- paste(
        unlist(lapply(seq_len(length(refs_split)), function(i) {
          refname <- refs_split_names[i]
          ref_doi <- refs_split_dois[i]
          .link_refs(refname = refname, ref_doi = ref_doi)
        })),
        collapse = " ")

      cell_html_refs <- paste0(
        tags$b("References: "),
        refs_content,
        tags$br(),
        tags$br()
      )


      # cellonts <- as.character(glossary_df[cell_type, "related_cellontology"])
      # cellonts_split <- unlist(strsplit(cellonts, split = ",", fixed = TRUE))
      # cellonts_content <- paste(
      #   unlist(lapply(cellonts_split, function(id) .link_cellontology(id))),
      #   collapse = " ")
      #
      # cell_html_cellontology <- paste0(
      #   tags$b("Cell ontology terms: "),
      #   cellonts_content,
      #   tags$br()
      # )

      ## TODO:
      # handle the part with including the umap plot, specifically for that subset
      if (display_plot) {
        img_location <- glossary_df[cell_type, "umap_png"]
        img_html <-
          paste0(
            tags$img(
              src = base64enc::dataURI(
                file = system.file("images", img_location, package = "HCATonsilData"),
                mime = "image/png"
              ),
              width = 500
            ),
            tags$hr()
          )
      } else {
        img_html <- ""
      }


      cell_html <- paste0(
        cell_html_celltype,
        img_html,
        cell_html_celldescription,
        cell_html_cellmarkers,
        cell_html_refs #,
        # cell_html_cellontology
      )




      if (output_to == "single_page") {
        # put this into a html file
        tmpfile <- tempfile()

        writeLines(cell_html, con = tmpfile)

        # and render this
        browseURL(rmarkdown::render(tmpfile))
        return(tmpfile)
      } else if (output_to == "html_to_embed") {
        return(cell_html)
      }


    } else {
      message("Cell type not found! Please select one of the following: ",
              paste0(slo_celltypes, collapse = "|"))
    }
  }
}
