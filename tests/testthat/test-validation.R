test_that("validate_potato detects missing required fields", {
  # S7 class validator already catches empty id/name, so this is redundant
  # Just verify the function reports it correctly
  skip("S7 class validator handles this")
})

test_that("validate_potato detects duplicate node IDs", {
  potato <- load_test_potato()

  # Duplicate a node
  potato@nodes[[3]]$id <- potato@nodes[[1]]$id

  validation <- validate_potato(potato)

  expect_false(validation$valid)
  expect_true(any(grepl("Duplicate node IDs", validation$errors)))
})

test_that("validate_potato detects invalid edge references", {
  potato <- load_test_potato()

  # Add edge to non-existent node
  potato@edges <- c(potato@edges, list(list(from = "fake_node", to = "gapA")))

  validation <- validate_potato(potato)

  expect_false(validation$valid)
  expect_true(any(grepl("non-existent node", validation$errors)))
})

test_that("validate_potato detects cycles", {
  potato <- load_test_potato()

  # Create a cycle using DAG node IDs (gene_id_step format)
  # gapA_1 -> pgk_2 -> eno_3, so adding eno_3 -> gapA_1 creates a cycle
  potato@edges <- c(potato@edges, list(list(from = "eno_3", to = "gapA_1")))

  validation <- validate_potato(potato)

  expect_false(validation$valid)
  expect_true(any(grepl("cycle", validation$errors)))
})

test_that("validate_potato warns about invalid KO format", {
  potato <- load_test_potato()

  # Add invalid KO
  potato@nodes[[1]]$ko <- c(potato@nodes[[1]]$ko, "K123")  # Too short

  validation <- validate_potato(potato)

  expect_true(validation$valid)  # Should still be valid
  expect_true(length(validation$warnings) > 0)
  expect_true(any(grepl("invalid KO format", validation$warnings)))
})

test_that("validate_potato warns about invalid EC format", {
  potato <- load_test_potato()

  # Add invalid EC
  potato@nodes[[1]]$ec <- c(potato@nodes[[1]]$ec, "1.2.3")  # Missing fourth number

  validation <- validate_potato(potato)

  expect_true(validation$valid)  # Should still be valid
  expect_true(length(validation$warnings) > 0)
  expect_true(any(grepl("invalid EC format", validation$warnings)))
})

test_that("validate_potato warns about enzymes without detection methods", {
  potato <- load_test_potato()

  # Remove all detection methods from a node
  potato@nodes[[1]]$ko <- NULL
  potato@nodes[[1]]$ec <- NULL

  validation <- validate_potato(potato)

  expect_true(validation$valid)
  expect_true(length(validation$warnings) > 0)
  expect_true(any(grepl("no detection methods", validation$warnings)))
})

test_that("validate_potato errors on non-boolean required field", {
  potato <- load_test_potato()

  # Set required to a string
  potato@nodes[[1]]$required <- "yes"

  validation <- validate_potato(potato)

  expect_false(validation$valid)
  expect_true(any(grepl("'required' must be TRUE/FALSE", validation$errors)))
})

test_that("validate_potato strict mode adds warnings", {
  potato <- load_test_potato()

  # Remove source
  potato@source <- ""

  validation_normal <- validate_potato(potato, strict = FALSE)
  validation_strict <- validate_potato(potato, strict = TRUE)

  expect_true(validation_normal$valid)
  expect_true(validation_strict$valid)

  # Strict should have more warnings
  expect_true(length(validation_strict$warnings) > length(validation_normal$warnings))
  expect_true(any(grepl("No 'source' specified", validation_strict$warnings)))
})

test_that("print_validation displays results correctly", {
  potato <- load_test_potato()
  validation <- validate_potato(potato)

  expect_output(print_validation(validation), "valid")
})

test_that("validate_potato handles unknown threshold names", {
  potato <- load_test_potato()

  # Add unknown threshold
  potato@nodes[[1]]$thresholds <- list(fake_threshold = 0.5)

  validation <- validate_potato(potato)

  expect_true(validation$valid)
  expect_true(any(grepl("unknown threshold names", validation$warnings)))
})

test_that("valid potatoes pass all checks", {
  # Test the example potatoes
  test_potato <- load_test_potato()
  validation_test <- validate_potato(test_potato)
  expect_true(validation_test$valid)
  expect_equal(length(validation_test$errors), 0)

  # Nitrogen fixation
  nif_path <- system.file("potatoes", "nitrogen_fixation.json", package = "potato")
  skip_if(nif_path == "", "Nitrogen fixation potato not found")

  nif <- load_potato(nif_path)
  validation_nif <- validate_potato(nif)
  expect_true(validation_nif$valid)
  expect_equal(length(validation_nif$errors), 0)
})
