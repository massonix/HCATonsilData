% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NBC_MBC_annotation_df}
\alias{NBC_MBC_annotation_df}
\title{Annotation dictionary for naive and memory B cells (NBC-MBC)}
\format{
A dataframe with 112478 NBC-MBC and 5 variables
\describe{
\item{barcode}{cell barcode}
\item{names_level_5_clusters_eta}{Cell annotation given by the B-cell annotation team in February 2022}
\item{annotation_20220619}{Cell annotation given by the B-cell annotation team June 2022 (preprint)}
\item{UMAP_1}{UMAP1 coordinates}
\item{UMAP_2}{UMAP2 coordinates}
}
}
\source{
\url{https://zenodo.org/record/8373756}
}
\description{
In version 1 of the atlas, we changed the clusters and annotation for NBC-MBC
prior to publishing the preprint. Since, there is not a 1:1 mapping between
clusters, here we provide the correspondence between the annotations in
February 2022 and the preprint (June 2022) for each NBC-MBC cell barcode.
}
\details{
This data is used by the updateAnnotation function in B cells.
}
