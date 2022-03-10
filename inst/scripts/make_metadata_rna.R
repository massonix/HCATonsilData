# This script generates the metadata_rna.csv file as requested in the following
# tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html


# Load packages
library(dplyr)
library(tidyr)
library(readr)
library(here)


# Read data
path_to_csv <- here("inst/scripts/tonsil_atlas_rds_files.csv")
files_df <- read.csv(file = path_to_csv, header = TRUE)
files_df <- files_df[files_df$dataset == "RNA", ]


# Generate metadata file for ExperimentHub package
files_df %>%
  mutate(outs = "coldata;counts;processed;rowdata;pca;harmony;umap") %>%
  separate_rows(outs, sep = ";") %>%
  mutate(
    outs2 = convs[outs],
    descs = convs2[outs],
    suffix = suffix[outs]) %>%
  View()
