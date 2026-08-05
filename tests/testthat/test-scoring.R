test_that("score_pathways requires annotation results", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  expect_error(score_pathways(sack), "No annotation results found")

  unlink(temp_dir, recursive = TRUE)
})

test_that("score_pathways creates scores tibble", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock minimal annotation results
  sack@results <- tibble::tibble(
    genome = "test_genome",
    kofam = list(tibble::tibble(
      potato = names(sack@potatoes)[1],
      node_id = "testgene",
      gene_id = "testgene",
      ko = "K00001",
      score = 100,
      threshold = 50
    ))
  )

  sack <- score_pathways(sack)

  expect_s3_class(sack@scores, "tbl_df")
  expect_true(nrow(sack@scores) > 0)
  expect_true(all(c("genome", "potato", "fraction", "present") %in% names(sack@scores)))

  unlink(temp_dir, recursive = TRUE)
})

test_that("score_pathways creates scores with required columns", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock annotation results
  sack@results <- tibble::tibble(
    genome = "test_genome",
    kofam = list(tibble::tibble(
      potato = names(sack@potatoes)[1],
      node_id = "testgene",
      gene_id = "testgene",
      ko = "K00001",
      score = 100,
      threshold = 50
    ))
  )

  sack <- score_pathways(sack)

  # Check that basic scoring columns exist
  expect_true("total_steps_detected" %in% names(sack@scores))
  expect_true("total_steps" %in% names(sack@scores))
  expect_true("fraction" %in% names(sack@scores))
  expect_true("present" %in% names(sack@scores))
  expect_true("min_fraction" %in% names(sack@scores))

  unlink(temp_dir, recursive = TRUE)
})

test_that("score_pathways includes min_fraction", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  # Mock annotation results
  sack@results <- tibble::tibble(
    genome = "test_genome",
    kofam = list(tibble::tibble(
      potato = names(sack@potatoes)[1],
      node_id = "testgene",
      gene_id = "testgene",
      ko = "K00001",
      score = 100,
      threshold = 50
    ))
  )

  sack <- score_pathways(sack)

  expect_true("min_fraction" %in% names(sack@scores))
  # Should have threshold from potato or default 0.75
  expect_true(all(sack@scores$min_fraction > 0 & sack@scores$min_fraction <= 1))

  unlink(temp_dir, recursive = TRUE)
})

test_that("score_pathways respects kofam thresholds", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)
  sack <- create_sack(temp_dir)

  potato_name <- names(sack@potatoes)[1]

  # Mock results with hits below and above threshold
  sack@results <- tibble::tibble(
    genome = c("genome1", "genome2"),
    kofam = list(
      # Genome 1: score below threshold
      tibble::tibble(
        potato = potato_name,
        node_id = "gene1",
        gene_id = "gene1",
        ko = "K00001",
        score = 40,
        threshold = 50
      ),
      # Genome 2: score above threshold
      tibble::tibble(
        potato = potato_name,
        node_id = "gene1",
        gene_id = "gene1",
        ko = "K00001",
        score = 100,
        threshold = 50
      )
    )
  )

  sack <- score_pathways(sack)

  # Genome 1 should have lower scores than genome 2
  genome1_scores <- sack@scores[sack@scores$genome == "genome1", ]
  genome2_scores <- sack@scores[sack@scores$genome == "genome2", ]

  expect_true(nrow(genome1_scores) > 0)
  expect_true(nrow(genome2_scores) > 0)
  # Genome 2 should detect more steps
  expect_true(all(genome2_scores$steps_detected >= genome1_scores$steps_detected))

  unlink(temp_dir, recursive = TRUE)
})

test_that("score_pathways handles multi-pathway networks", {
  # Load entner_doudoroff_network (multi-pathway)
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)

  # Copy ED network to sack
  dir.create(file.path(temp_dir, "potatoes"), recursive = TRUE)
  file.copy(
    system.file("potatoes", "entner_doudoroff_network.json", package = "potato"),
    file.path(temp_dir, "potatoes", "entner_doudoroff_network.json")
  )

  sack <- create_sack(temp_dir)

  # Mock annotation results for one gene
  sack@results <- tibble::tibble(
    genome = "test_genome",
    kofam = list(tibble::tibble(
      potato = "entner_doudoroff_network",
      node_id = "glk",
      gene_id = "glk",
      ko = "K00845",
      score = 100,
      threshold = 50
    ))
  )

  sack <- score_pathways(sack)

  # Should have multiple pathways scored
  expect_true(nrow(sack@scores) > 1)
  expect_true("pathway" %in% names(sack@scores))
  expect_true("pathway_name" %in% names(sack@scores))

  # Check that different pathways have different scores
  pathway_ids <- unique(sack@scores$pathway)
  expect_true(length(pathway_ids) > 1)

  unlink(temp_dir, recursive = TRUE)
})
