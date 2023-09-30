# This script generates the metadata.csv file for the Spatial data as requested in the following
# tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html
# Note that this script is inspired in https://github.com/fmicompbio/TabulaMurisSenisData/blob/master/inst/scripts/make-metadata-tabula-muris-senis-droplet.R

# Load packages
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(here)
library(glue)
library(purrr)


# Get data slots
data_dir <- here("inst/scripts/HCATonsilData/2.0/Spatial")
slots <- list.files(data_dir)


# Generate metadata file for ExperimentHub package
slots_no_image <- slots %>%
  str_subset("(image|scale)", negate = TRUE)
slots_names <- slots_no_image %>%
  str_split("_") %>%
  map_chr(\(x) str_subset(x, "(.rds|.h5)")) %>%
  str_remove("(.rds|.h5)")
titles <- c(counts = "counts", logcounts = "processed counts", coldata = "colData",
            rowdata = "rowData", pca = "PCA", harmony = "HARMONY", umap = "UMAP")
descs <- c(counts = "Count matrix", logcounts = "Processed count matrix", coldata = "Cell metadata",
           rowdata = "Gene annotation", pca = "PCA representation", harmony = "HARMONY representation",
           umap = "UMAP representation")
rdclass <- c(counts = "H5File", logcounts = "H5File", coldata = "DFrame",
             rowdata = "DFrame", pca = "matrix", harmony = "matrix",
             umap = "matrix")
names(slots_no_image) <- slots_names
metadata_df1 <- data.frame(
  Title = glue("Tonsil Atlas Visium data {titles[slots_names]}"),
  Description = glue("{descs[slots_names]} for the Tonsil Cell Atlas Visium dataset"),
  RDataPath = glue("HCATonsilData/2.0/Spatial/{slots_no_image[slots_names]}"),
  RDataClass = rdclass[slots_names],
  DispatchClass = ifelse(rdclass[slots_names] == "H5File", "H5File", "Rds")
)


# Handle images and scales
slots_image <- str_subset(slots, "(image|scale)")
sample_ids <- slots_image %>%
  str_remove("(Spatial_image_|Spatial_scale_)") %>%
  str_remove(".rds")
file_type <- str_extract(slots_image, "(scale|image)")
file_type[file_type == "image"] <- "H&E staining"
metadata_df2 <- data.frame(
  Title = glue("Tonsil Atlas Visium {file_type} ({sample_ids})"),
  Description = glue("{file_type} for sample {sample_ids}"),
  RDataPath = glue("HCATonsilData/2.0/Spatial/{slots_image}"),
  RDataClass = ifelse(file_type == "scale", "scalefactors", "raster"),
  DispatchClass = "Rds"
)
metadata_df <- bind_rows(metadata_df1, metadata_df2)
metadata_df <- metadata_df %>%
  mutate(
    BiocVersion = "3.18",
    Genome = "GRCh38",
    SourceType = "HDF5",
    SourceUrl = "https://zenodo.org/record/8373756",
    SourceVersion = "2.0",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = NA,
    DataProvider = "BCLL@las",
    Maintainer = "Ramon Massoni-Badosa <ramonmassoni@gmail.com>"
) %>%
  select(Title, Description, RDataPath, BiocVersion, Genome,
         SourceType, SourceUrl, SourceVersion, Species, TaxonomyId,
         Coordinate_1_based, DataProvider, Maintainer, RDataClass,
         DispatchClass)
rownames(metadata_df) <- NULL


# Write
write_delim(
  metadata_df,
  file = here("inst/extdata/metadata-tonsil-atlas-spatial_v2.csv"),
  delim = ","
)
