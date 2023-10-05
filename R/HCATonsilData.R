#' Access the Tonsil Atlas data (RNA, ATAC, Multiome, CITE, Spatial)
#'
#' The data was downloaded from Zenodo
#' https://zenodo.org/record/8373756
#'
#' @param assayType One of 'RNA', 'ATAC', 'Multiome', 'CITE' or 'Spatial'.
#' @param cellType A character vector of length 1 with the desired cell type.
#'    A list of available cell types can be obtained using
#'    \code{listCellTypes(assay_type)}.
#' @param processedCounts Logical scalar. If \code{TRUE}, include the processed
#'   (normalized) counts in addition to the raw counts in the
#'   SingleCellExperiment object.
#' @param version Version of the tonsil atlas data to retrieve: "1.0" (preprint)
#'   or "2.0" (publication, default)
#'
#' @return A \linkS4class{SingleCellExperiment} object for the cellType
#'   requested. For scATAC-seq, Multiome and CITE we provide the instructions
#'   for downloading the Seurat objects in Zenodo (see vignette)
#'
#' @author Ramon Massoni-Badosa

#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom HDF5Array HDF5Array
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim reducedDim<-
#'   reducedDims
#' @importFrom SpatialExperiment SpatialExperiment SpatialImage
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @examples
#' # retrieve the epithelial scRNA-seq dataset
#' sce_epithelial <- HCATonsilData(
#'   assayType = "RNA",
#'   cellType = "epithelial"
#' )
#' sce_epithelial

HCATonsilData <- function(assayType = c("RNA", "ATAC", "CITE", "Spatial"),
                          cellType = "All",
                          version = "2.0",
                          processedCounts = TRUE) {

  # Check validity of input arguments
  assayType <- match.arg(assayType)
  version <- match.arg(version, c("1.0", "2.0"))
  cellTypes <- listCellTypes(assayType, version)
  cellType <- match.arg(cellType, cellTypes)

  # Point users to the vignette to download scATAC-seq, Multiome or CITE-seq data
  if (assayType %in% c("ATAC", "Multiome", "CITE")) {
    message(
      sprintf(
        "Check the corresponding section in the vignette to download %s data",
        assayType
      )
    )
    message("Please run browseVignettes('HCATonsilData')")
  }

  # Download scRNA-seq data and return a SingleCellExperiment if assay is RNA
  if (assayType == "RNA") {
    # Initialize ExperimentHub object
    eh <- ExperimentHub::ExperimentHub()
    host <- file.path("HCATonsilData", version, assayType)

    # Define paths to files in ExperimentHub
    suffixes <- c("counts.h5", "processed.h5", "rowdata.rds", "coldata.rds",
                  "pca.rds", "harmony.rds", "umap.rds")
    names(suffixes) <- c("counts", "processedCounts", "rowData", "colData",
                         "pca", "harmony", "umap")
    filePaths <- file.path(host, paste(cellType, assayType, suffixes, sep = "_"))
    names(filePaths) <- names(suffixes)
    for (x in filePaths) {
      if (sum(x == eh$rdatapath) > 1) {
        stop(x, " matched more than one entry!")
      }
      if (sum(x == eh$rdatapath) == 0) {
        stop(x, " is not present in ExperimentHub!")
      }
    }

    # Create SingleCellExperiment
    cts <- HDF5Array::HDF5Array(
      eh[eh$rdatapath == filePaths["counts"]][[1]],
      name = "counts"
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = cts),
      rowData = eh[eh$rdatapath == filePaths["rowData"]][[1]],
      colData = eh[eh$rdatapath == filePaths["colData"]][[1]]
    )
    SingleCellExperiment::reducedDims(sce) <- list(
      PCA = eh[eh$rdatapath == filePaths["pca"]][[1]],
      UMAP = eh[eh$rdatapath == filePaths["umap"]][[1]]
    )
    if (sum(filePaths["harmony"] == eh$rdatapath) == 1) {
      harmony <- eh[eh$rdatapath == filePaths["harmony"]][[1]]
      SingleCellExperiment::reducedDims(sce)[["HARMONY"]] <- harmony
    }
    if (processedCounts) {
      prccts <- HDF5Array::HDF5Array(
        eh[eh$rdatapath == filePaths["processedCounts"]][[1]],
        name = "processed"
      )
      SummarizedExperiment::assay(sce, "logcounts", withDimnames = FALSE) <- prccts
    }
    
    if (version == "1.0") {
      # Update annotations
      sce <- updateAnnotation(sce = sce, "20220215", "20220619")
      sce <- updateAnnotation(sce = sce, "20220619", "20230508")
  
      # For NBC-MBC, we reclustered them from 20220215 to 20220619, which means
      # that there is not a 1:1 mapping between annotations. Let us reannotate
      # using a pre-loaded dataframe. Also, let's change the UMAP coords
      # to match the paper
      if (cellType == "NBC-MBC") {
        sce <- sce[, NBC_MBC_annotation_df$barcode]
        sce$annotation_20220619 <- HCATonsilData::NBC_MBC_annotation_df$annotation_20220619
        umap_df <- as.matrix(
          HCATonsilData::NBC_MBC_annotation_df[, c("UMAP_1", "UMAP_2")]
        )
        reducedDim(sce, "UMAP") <- umap_df
      }
  
      # Similarly, let us update NBC/MBC annotation in case cellType = All
      if (cellType == "All") {
        annot <- sce$annotation_20220619
        names(annot) <- colnames(sce)
        annot[HCATonsilData::NBC_MBC_annotation_df$barcode] <- HCATonsilData::NBC_MBC_annotation_df$annotation_20220619
        sce$annotation_20220619 <- annot
      }
    }
  }
    
  # Get Visium data and return SpatialExperiment if assaytype is Spatial
  if (assayType == "Spatial") {
    if (version == "1.0") {
      stop("Spatial transcriptomics data is only available for version 2.0")
    }
    
    # Initialize ExperimentHub object
    eh <- ExperimentHub::ExperimentHub()
    host <- file.path("HCATonsilData", version, assayType)
    
    # Get raw and processed counts
    cts <- HDF5Array::HDF5Array(
      eh[eh$rdatapath == file.path(host, "Spatial_assay_counts.h5")][[1]],
      name = "counts"
    )
    as <- list(counts = cts)
    if (processedCounts) {
      prccts <- HDF5Array::HDF5Array(
        eh[eh$rdatapath == file.path(host, "Spatial_assay_logcounts.h5")][[1]],
        name = "logcounts"
      )
      as[["logcounts"]] <- prccts
    }
    
    # Get dimensionality reductions
    dr_filt <- grepl(file.path(host, "Spatial_dimred"), eh$rdatapath)
    dr_paths <- eh$rdatapath[dr_filt]
    names(dr_paths) <- sub(".*/Spatial_dimred_(\\w+)\\.rds", "\\1", dr_paths)
    dr <- lapply(dr_paths, \(x) eh[eh$rdatapath == x][[1]])
    names(dr) <- toupper(names(dr))
    
    # Get cell and gene metadata
    rowdata <- eh[eh$rdatapath == file.path(host, "Spatial_rowdata.rds")][[1]]
    coldata <- eh[eh$rdatapath == file.path(host, "Spatial_coldata.rds")][[1]]
    
    # Get images
    img_filt <- grepl(file.path(host, "Spatial_image"), eh$rdatapath)
    sfs_filt <- grepl(file.path(host, "Spatial_scale"), eh$rdatapath)
    img_paths <- eh$rdatapath[img_filt]
    sfs_paths <- eh$rdatapath[sfs_filt]
    names(img_paths) <- sub(".*/Spatial_image_(\\w+)\\.rds", "\\1", img_paths)
    names(sfs_paths) <- sub(".*/Spatial_scale_(\\w+)\\.rds", "\\1", sfs_paths)
    sfs_paths <- sfs_paths[names(img_paths)] # ensure images and scale are ordered equally
    images <- mapply(\(img_path, sfs_path, nm) {
      img <- eh[eh$rdatapath == img_path][[1]]
      sfs <- eh[eh$rdatapath == sfs_path][[1]]
      df <- DataFrame(
        sample_id = nm,
        image_id = "lowres",
        data = I(list(SpatialImage(img))),
        scaleFactor = sfs$lowres
      )
      df
    }, img_path=img_paths, sfs_path = sfs_paths, nm=names(img_paths), SIMPLIFY = TRUE)
    images <- do.call(rbind, images)

    # Create SpatialExperiment object
    sce <- SpatialExperiment(
      assays = as,
      mainExpName = "Spatial",
      reducedDims = dr,
      imgData = images,
      rowData = rowdata,
      colData = coldata
    )
  }

  # Return
  return(sce)
}


#' #' Information on the Tonsil Data
#' #'
#' #' Get information on the individual datasets of the Tonsil Atlas data
#' #'
#' #' @inheritParams HCATonsilData
#' #'
#' #' @return An info message is generated, and its content is returned invisibly.
#' #'
#' #' @author Federico Marini
#' 
#' #' @export
#' #'
#' #' @examples
#' #' HCATonsilDataInfo(assayType = "RNA", cellType = "epithelial")
#' HCATonsilDataInfo <- function(assayType = c("RNA", "ATAC", "CITE", "Spatial"),
#'                               cellType = listCellTypes(assayType = assayType)) {
#' 
#'   # Check validity of input arguments
#'   assayType <- match.arg(assayType)
#'   cellType <- match.arg(cellType)
#' 
#'   # Initialize ExperimentHub object
#'   eh <- ExperimentHub::ExperimentHub()
#'   host <- file.path("HCATonsilData", "1.0", assayType)
#' 
#'   # Define paths
#'   suffixes <- c("counts.h5", "processed.h5", "rowdata.rds", "coldata.rds",
#'                 "pca.rds", "harmony.rds", "umap.rds")
#'   names(suffixes) <- c("counts", "processedCounts", "rowData", "colData",
#'                        "pca", "harmony", "umap")
#'   filePaths <- vapply(suffixes, \(.) {
#'     x <- file.path(host, paste(cellType, assayType, ., sep = "_"))
#'     x
#'   }, character(1))
#'   for (x in filePaths) {
#'     if (sum(x == eh$rdatapath) > 1) {
#'       stop("Input matched more than one entry!")
#'     }
#'   }
#' 
#'   # Fetch info instead of creating full sce object
#'   rd_sce <- eh[eh$rdatapath == filePaths["rowData"]][[1]]
#'   cd_sce <- eh[eh$rdatapath == filePaths["colData"]][[1]]
#' 
#'   # buildup the information string
#'   info_sce <- sprintf(
#'     "Dataset: %s (assay selected: %s)\nMeasurements available for %d genes in %d cells",
#'     cellType,
#'     assayType,
#'     nrow(rd_sce),
#'     nrow(cd_sce)
#'   )
#' 
#'   message(info_sce)
#' 
#'   return(invisible(info_sce))
#' }
