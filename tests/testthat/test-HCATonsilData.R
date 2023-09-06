test_that("downloading RNA data works", {
  # Download data
  sce <- HCATonsilData(
    assayType = "RNA",
    cellType = "epithelial",
    processedCounts = TRUE
  )

  # Test slots and classes
  expect_s4_class(sce, "SingleCellExperiment")
  expect_s4_class(counts(sce), "DelayedMatrix")
  expect_s4_class(SingleCellExperiment::logcounts(sce), "DelayedMatrix")
  expect_equal(length(SingleCellExperiment::reducedDims(sce)), 3)
  expect_equal(names(SingleCellExperiment::reducedDims(sce)), c("PCA", "UMAP", "HARMONY"))
  expect_s4_class(SingleCellExperiment::colData(sce), "DFrame")
  expect_s4_class(SingleCellExperiment::rowData(sce), "DFrame")

  # Test dimensions
  expect_equal(dim(sce), c(37378, 277))
  expect_equal(dim(counts(sce)), c(37378, 277))
  expect_equal(dim(SingleCellExperiment::logcounts(sce)), c(37378, 277))
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, "PCA")), c(277, 50))
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, "UMAP")), c(277, 2))
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, "HARMONY")), c(277, 50))
  expect_equal(dim(SingleCellExperiment::rowData(sce)), c(37378, 3))
  expect_equal(dim(SingleCellExperiment::colData(sce)), c(277, 31))

  # Test variables colData()
  vars <- c("barcode", "donor_id", "gem_id", "library_name", "assay", "sex",
            "age", "age_group", "hospital", "UMAP_1_level_1", "UMAP_2_level_1",
            "nCount_RNA", "nFeature_RNA", "pct_mt", "pct_ribosomal",
            "is_hashed", "pDNN_hashing", "pDNN_scrublet", "pDNN_union",
            "scrublet_doublet_scores", "scrublet_predicted_doublet", "S.Score",
            "G2M.Score", "Phase", "CC.Difference", "annotation_level_1",
            "annotation_figure_1", "annotation_20220215", "UMAP_1_20220215",
            "UMAP_2_20220215")
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
  cellTypes <- HCATonsilData::listCellTypes("RNA")
  out <- all(sapply(cellTypes, \(.) any(grepl(., eh$rdatapath))))
  expect_true(out)
})

test_that("Info retrieval", {
  expect_message(
    info_sce <- HCATonsilDataInfo(assayType = "RNA", cellType = "epithelial")
  )
  expect_type(info_sce, "character")
  expect_true(grepl(pattern = "277 cells", info_sce))
})


test_that("Updating annotations works", {
  sce <- HCATonsilData(
    assayType = "RNA",
    cellType = "CD8-T",
    processedCounts = FALSE
  )
  sce <- updateAnnotation(
    sce = sce,
    refAnnotation = "20220215",
    newAnnotation = "20220619"
  )
  oldAnnot <- sce$annotation_20220215
  newAnnot <- sce$annotation_20220619
  expect_true(
    all(colnames(sce)[oldAnnot == "CXCR6+ RM CD8 T"] ==
        colnames(sce)[newAnnot == "RM CD8 activated T"])
  )
})
