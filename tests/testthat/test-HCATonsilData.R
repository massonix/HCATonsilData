data("annotations_dictionary")
test_that("downloading RNA data works for both versions", {
  # Version 1
  sce_v1 <- HCATonsilData(
    assayType = "RNA",
    cellType = "epithelial",
    processedCounts = TRUE,
    version = "1.0"
  )
  
  # Version 2
  sce <- HCATonsilData(
    assayType = "RNA",
    cellType = "epithelial",
    processedCounts = TRUE,
    version = "2.0"
  )

  # Test slots and classes
  expect_s4_class(sce_v1, "SingleCellExperiment")
  expect_s4_class(sce, "SingleCellExperiment")
  expect_s4_class(counts(sce), "DelayedMatrix")
  expect_s4_class(SingleCellExperiment::logcounts(sce), "DelayedMatrix")
  expect_equal(length(SingleCellExperiment::reducedDims(sce)), 3)
  expect_equal(names(SingleCellExperiment::reducedDims(sce)), c("PCA", "UMAP", "HARMONY"))
  expect_s4_class(SingleCellExperiment::colData(sce), "DFrame")
  expect_s4_class(SingleCellExperiment::rowData(sce), "DFrame")

  # Test dimensions
  expect_equal(dim(sce), c(37378, 396))
  expect_equal(dim(counts(sce)), c(37378, 396))
  expect_equal(dim(SingleCellExperiment::logcounts(sce)), c(37378, 396))
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, "PCA")), c(396, 50))
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, "UMAP")), c(396, 2))
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, "HARMONY")), c(396, 50))
  expect_equal(dim(SingleCellExperiment::rowData(sce)), c(37378, 3))
  expect_equal(dim(SingleCellExperiment::colData(sce)), c(396, 39))

  # Test variables colData()
  vars <- c("barcode", "donor_id", "gem_id", "library_name", "assay", "sex",
            "age","age_group", "hospital", "cohort_type", "cause_for_tonsillectomy",
            "is_hashed", "preservation", "nCount_RNA", "nFeature_RNA", "pct_mt",
            "pct_ribosomal", "pDNN_hashing", "pDNN_scrublet", "pDNN_union",
            "scrublet_doublet_scores", "S.Score", "G2M.Score", "Phase",
            "scrublet_predicted_doublet", "doublet_score_scDblFinder",
            "annotation_level_1", "annotation_level_1_probability",
            "annotation_figure_1", "annotation_20220215", "annotation_20220619",
            "annotation_20230508", "annotation_20230508_probability",
            "UMAP_1_level_1", "UMAP_2_level_1", "UMAP_1_20220215", "UMAP_2_20220215",
            "UMAP_1_20230508", "UMAP_2_20230508")
  expect_true(all(vars %in% colnames(SingleCellExperiment::colData(sce))))

  # Test that HCATonsilData generates errors if inputs are incorrect
  expect_error(HCATonsilData(assayType = "DNA"), regexp = NULL)
  expect_error(
    HCATonsilData(assayType = "RNA", cellType = "neuron"),
    regexp = NULL
  )

  # Test that logcounts are not included if processedCounts = FALSE
  sce2 <- HCATonsilData(
    assayType = "RNA",
    cellType = "epithelial",
    processedCounts = FALSE
  )
  expect_error(SingleCellExperiment::logcounts(sce2), regexp = NULL)

  # Test that

})


test_that("All cellTypes exist in ExperimentHub", {
  eh <- ExperimentHub::ExperimentHub()
  cellTypes <- HCATonsilData::listCellTypes("RNA", version = "2.0")
  out <- all(vapply(cellTypes, \(.) any(grepl(., eh$rdatapath)), logical(1)))
  expect_true(out)
  cellTypes <- HCATonsilData::listCellTypes("RNA", version = "1.0")
  out <- all(vapply(cellTypes, \(.) any(grepl(., eh$rdatapath)), logical(1)))
  expect_true(out)
})

# test_that("Info retrieval", {
#   expect_message(
#     info_sce <- HCATonsilDataInfo(assayType = "RNA", cellType = "epithelial")
#   )
#   expect_type(info_sce, "character")
#   expect_true(grepl(pattern = "277 cells", info_sce))
# })


test_that("Updating annotations works", {
  sce <- HCATonsilData(
    assayType = "RNA",
    cellType = "myeloid",
    processedCounts = FALSE,
    version = "1.0"
  )
  # 20220215 --> 20220619
  sce <- updateAnnotation(
    sce = sce,
    refAnnotation = "20220215",
    newAnnotation = "20220619"
  )
  oldAnnot <- sce$annotation_20220215
  newAnnot <- sce$annotation_20220619
  expect_true(
    all(colnames(sce)[oldAnnot == "SELENOP FUCA1 PTGDS macrophages"] ==
        colnames(sce)[newAnnot == "SELENOP Slancytes"])
  )
  
  # 20220619 --> 20230508
  sce <- updateAnnotation(
    sce = sce,
    refAnnotation = "20220619",
    newAnnotation = "20230508"
  )
  oldAnnot <- sce$annotation_20220619
  newAnnot <- sce$annotation_20230508
  expect_true(
    all(colnames(sce)[oldAnnot == "SELENOP FUCA1 PTGDS macrophages"] ==
          colnames(sce)[newAnnot == "SELENOP Slan-like"])
  )
})
