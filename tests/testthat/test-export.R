test_that("get_gene_results returns empty tibble for empty sack", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)
  sack <- create_sack(temp_dir)

  result <- get_gene_results(sack)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)

  unlink(temp_dir, recursive = TRUE)
})

test_that("get_gene_results includes passed column", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock some kofam results
  sack@results <- tibble::tibble(
    genome = c("test_genome"),
    kofam = list(tibble::tibble(
      potato = "test_potato",
      node_id = "geneA",
      gene_id = "geneA",
      ko = "K00001",
      score = 100,
      threshold = 50
    ))
  )

  result <- get_gene_results(sack)

  expect_true("passed" %in% names(result))
  expect_true("tool" %in% names(result))
  expect_equal(result$passed[1], TRUE)  # score >= threshold
  expect_equal(result$tool[1], "kofam")

  unlink(temp_dir, recursive = TRUE)
})

test_that("get_gene_results handles multiple tools", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock results from multiple tools
  sack@results <- tibble::tibble(
    genome = c("test_genome"),
    kofam = list(tibble::tibble(
      potato = "test_potato",
      node_id = "geneA",
      gene_id = "geneA",
      ko = "K00001",
      score = 100,
      threshold = 50
    )),
    blast = list(tibble::tibble(
      potato = "test_potato",
      node_id = "geneB",
      gene_id = "geneB",
      blast_id = "ref1",
      evalue = 1e-20,
      bitscore = 200
    )),
    hmm = list(tibble::tibble(
      potato = "test_potato",
      node_id = "geneC",
      gene_id = "geneC",
      profile = "PF00001",
      score = 50,
      evalue = 1e-15,
      tc_threshold = 40
    ))
  )

  result <- get_gene_results(sack)

  expect_equal(nrow(result), 3)
  expect_true(all(c("kofam", "blast", "hmm") %in% result$tool))
  expect_true(all(result$passed))  # All should pass thresholds

  unlink(temp_dir, recursive = TRUE)
})

test_that("get_pathway_scores returns empty tibble for empty sack", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)
  sack <- create_sack(temp_dir)

  result <- get_pathway_scores(sack)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)

  unlink(temp_dir, recursive = TRUE)
})

test_that("get_pathway_scores includes potato_hash", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock some scores (for single-pathway potato)
  potato_name <- names(sack@potatoes)[1]
  sack@scores <- tibble::tibble(
    genome = "test_genome",
    potato = potato_name,
    potato_name = "Test Potato",
    steps_detected = 5,
    steps_total = 8,
    fraction = 0.625,
    min_fraction = 0.75,
    present = FALSE,
    steps_detected_essential = 5,
    steps_total_essential = 5,
    fraction_essential = 1.0,
    present_essential = TRUE
  )

  result <- get_pathway_scores(sack)

  expect_true("potato_hash" %in% names(result))
  expect_false("potato" %in% names(result))  # Should be removed

  unlink(temp_dir, recursive = TRUE)
})

test_that("get_pathway_scores includes essential metrics", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock scores with essential metrics
  sack@scores <- tibble::tibble(
    genome = "test_genome",
    potato = "test_potato",
    potato_name = "Test Potato",
    steps_detected = 5,
    steps_total = 10,
    fraction = 0.5,
    min_fraction = 0.75,
    present = FALSE,
    steps_detected_essential = 5,
    steps_total_essential = 5,
    fraction_essential = 1.0,
    present_essential = TRUE
  )

  result <- get_pathway_scores(sack)

  expect_true("steps_detected_essential" %in% names(result))
  expect_true("fraction_essential" %in% names(result))
  expect_true("present_essential" %in% names(result))
  expect_equal(result$fraction_essential[1], 1.0)

  unlink(temp_dir, recursive = TRUE)
})

test_that("get_pathway_scores includes min_fraction threshold", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock scores
  sack@scores <- tibble::tibble(
    genome = "test_genome",
    potato = "test_potato",
    potato_name = "Test Potato",
    steps_detected = 5,
    steps_total = 10,
    fraction = 0.5,
    min_fraction = 0.8,
    present = FALSE,
    steps_detected_essential = NA_integer_,
    steps_total_essential = NA_integer_,
    fraction_essential = NA_real_,
    present_essential = NA
  )

  result <- get_pathway_scores(sack)

  expect_true("min_fraction" %in% names(result))
  expect_equal(result$min_fraction[1], 0.8)

  unlink(temp_dir, recursive = TRUE)
})
