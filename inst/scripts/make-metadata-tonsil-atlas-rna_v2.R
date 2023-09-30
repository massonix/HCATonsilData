# This script generates the metadata_rna.csv file as requested in the following
# tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html
# Note that this script is inspired in https://github.com/fmicompbio/TabulaMurisSenisData/blob/master/inst/scripts/make-metadata-tabula-muris-senis-droplet.R

# Load packages
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(here)
library(glue)


# Read data
data_dir <- here("inst/scripts/HCATonsilData/2.0/RNA")
slots <- list.files(data_dir)
cell_types <- slots %>%
  str_split("_") %>%
  map_chr(1) %>%
  unique()
df <- data.frame(cell_type = cell_types)


# Generate metadata file for ExperimentHub package
convs <- c(coldata = "colData", counts = "counts", processed = "processed counts",
           rowdata = "rowData", pca = "PCA", harmony = "HARMONY", umap = "UMAP")
convs2 <- c(coldata = "Cell metadata", counts = "Count matrix",
            processed = "Processed count matrix", rowdata = "Gene annotation",
            pca = "PCA representation", harmony = "HARMONY representation",
            umap = "UMAP representation")
suffix <- c(coldata = ".rds", counts = ".h5", processed = ".h5", rowdata = ".rds",
            pca = ".rds", harmony = ".rds", umap = ".rds")
rdclass <- c(coldata = "DFrame", counts = "H5File", processed = "H5File",
             rowdata = "DFrame", pca = "matrix", harmony = "matrix", umap = "matrix")
out_df <- df %>%
  mutate(outs = "coldata;counts;processed;rowdata;pca;harmony;umap") %>%
  separate_rows(outs, sep = ";") %>%
  mutate(
    outs2 = convs[outs],
    descs = convs2[outs],
    suffix = suffix[outs]
  ) %>%
  mutate(cell_type = str_replace(cell_type, "(/| )", "_")) %>%
  mutate(
    Title = glue("Tonsil Atlas RNA {cell_type} {outs2}"),
    Description = glue("{descs} for the Tonsil Cell Atlas {cell_type} scRNA-seq/Multiome dataset"),
    RDataPath = glue("HCATonsilData/2.0/RNA/{cell_type}_{dataset}_{outs}{suffix}"),
    BiocVersion = "3.18",
    Genome = "GRCh38",
    SourceType = "HDF5",
    SourceUrl = "https://zenodo.org/record/8373756",
    SourceVersion = "2.0",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = NA,
    DataProvider = "BCLL@las",
    Maintainer = "Ramon Massoni-Badosa <ramonmassoni@gmail.com>",
    RDataClass = rdclass[outs],
    DispatchClass = ifelse(suffix == ".h5", "H5File", "Rds")
  ) %>%
  select(Title, Description, RDataPath, BiocVersion, Genome,
         SourceType, SourceUrl, SourceVersion, Species, TaxonomyId,
         Coordinate_1_based, DataProvider, Maintainer, RDataClass,
         DispatchClass)

# Check that all files exist
if (all(str_remove(out_df$RDataPath, "HCATonsilData/2.0/RNA/") %in% slots)) {
  print("All files exist and are ready to be uploaded")
} else {
  print("Missing files!")
}


# Write
write_delim(
  out_df,
  file = here("inst/extdata/metadata-tonsil-atlas-rnal-v2.csv"),
  delim = ","
)
