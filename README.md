# HCATonsilData

HCATonsilData is an R/ExperimentHub package that provides easy access to single-cell RNA-seq (scRNA-seq), single-cell ATAC-seq (scATAC-seq), 10X Multiome, CITE-seq and spatial transcriptomics data (Visium) derived from the tonsil cell atlas project. It is inspired in the [TabulaMurisSenisData](https://github.com/fmicompbio/TabulaMurisSenisData/blob/master/README.md) package. For now, it can be installed as follows

``` {r}
devtools::install_github("massonix/HCATonsilData", build_vignettes = TRUE)
```

This will compile a vignette which documents the whole dataset:

``` {r}
browseVignettes("HCATonsilData")
```

# Requirements

The data was uploaded to ExperimentHub in the release 3.15 of Bioconductor, Thus, you need to have
BioC >= 3.15 to be able to use the package. In addition, you will need R >=4.2.


# Available data

## RNA


## ATAC


## CITE


## Spatial
