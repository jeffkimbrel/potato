test_that("summarize_annotation returns correct structure", {
  skip_if_not_installed("ggplot2")

  # This would require a real annotated sack, so skip for now
  skip("Requires fully annotated sack - integration test needed")

  # When implemented:
  # result <- summarize_annotation(sack)
  # expect_type(result, "list")
  # expect_named(result, c("status", "summary", "plot", "messages"))
  # expect_type(result$status, "list")
  # expect_s3_class(result$summary, "data.frame")
  # expect_s3_class(result$messages, "data.frame")
})

test_that("summarize_scoring returns correct structure", {
  skip_if_not_installed("ggplot2")

  # This would require a real scored sack, so skip for now
  skip("Requires fully scored sack - integration test needed")

  # When implemented:
  # result <- summarize_scoring(sack)
  # expect_type(result, "list")
  # expect_named(result, c("status", "summary", "plot", "messages"))
  # expect_type(result$status, "list")
  # expect_s3_class(result$summary, "data.frame")
  # expect_s3_class(result$messages, "data.frame")
})

test_that("get_annotation_details filters correctly", {
  skip("Requires fully annotated sack - integration test needed")

  # When implemented:
  # details <- get_annotation_details(sack)
  # expect_s3_class(details, "data.frame")
  # expect_true("database" %in% names(details))
  # expect_true("threshold_message" %in% names(details))

  # Test filtering by tool
  # kofam_only <- get_annotation_details(sack, tool = "kofam")
  # expect_true(all(kofam_only$database %in% kofam_db_names))
})

test_that("get_scoring_details filters correctly", {
  skip("Requires fully scored sack - integration test needed")

  # When implemented:
  # details <- get_scoring_details(sack)
  # expect_s3_class(details, "data.frame")
  # expect_true("present_via_marker" %in% names(details))
  # expect_true("marker_genes_found" %in% names(details))

  # Test present_only filter
  # present_only <- get_scoring_details(sack, present_only = TRUE)
  # expect_true(all(present_only$present))
})

test_that("get_detected_genes returns correct format", {
  skip("Requires fully scored sack - integration test needed")

  # When implemented:
  # genes_vec <- get_detected_genes(sack, genome = "genome001", potato = "test_glycolysis", format = "vector")
  # expect_type(genes_vec, "character")

  # genes_df <- get_detected_genes(sack, genome = "genome001", potato = "test_glycolysis", format = "dataframe")
  # expect_s3_class(genes_df, "data.frame")
  # expect_true("gene_id" %in% names(genes_df))
  # expect_true("is_marker" %in% names(genes_df))
})
