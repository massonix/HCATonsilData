test_that("loading glossary", {
  glossary_df <- TonsilData_glossary()
  expect_s3_class(glossary_df, "data.frame")

  expect_true(ncol(glossary_df) == 7)
})


test_that("Picking up info on cell type", {

  expect_message(TonsilData_cellinfo("preB"))
  expect_message(TonsilData_cellinfo())
  expect_message(TonsilData_cellinfo("Some cell not in there"))

  expect_true(
    is(TonsilData_cellinfo_html("Mature IgG+ PC"), "character")
  )

  expect_true(
    is(TonsilData_cellinfo_html("Mature IgG+ PC", output_to = "html_to_embed"), "character")
  )

  expect_message(
    TonsilData_cellinfo_html()
  )

  expect_message(
    TonsilData_cellinfo_html("A weird cell type")
  )


})

