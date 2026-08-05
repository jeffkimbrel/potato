test_that("load_potato handles multi-pathway networks", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  expect_true(S7::S7_inherits(potato, Potato))
  expect_true(!is.null(potato@edges))
  expect_true(is.list(potato@edges))
  expect_true(length(potato@edges) > 1)  # Multiple pathways

  # Check pathway structure
  pathway <- potato@edges[[1]]
  expect_true(!is.null(pathway$name))
  expect_true(!is.null(pathway$type))
  expect_true(!is.null(pathway$nodes))
  expect_true(!is.null(pathway$edges))
})

test_that("validate_potato works with multi-pathway networks", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  validation <- validate_potato(potato)

  expect_true(validation$valid)
  expect_equal(length(validation$errors), 0)
})

test_that("print_potato works with multi-pathway networks", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  # Should not error when printing
  expect_no_error(print_potato(potato))
})

test_that("build_potato_graph works with multi-pathway networks", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  g <- build_potato_graph(potato)

  expect_s3_class(g, "igraph")
  expect_true(igraph::vcount(g) > 0)
  expect_true(igraph::ecount(g) > 0)

  # Node names should be gene IDs (not gene_step format)
  node_names <- igraph::V(g)$name
  expect_false(any(grepl("_\\d+$", node_names)))
})

test_that("multi-pathway networks have shared genes", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  # Get genes from each pathway
  pathway_genes <- lapply(potato@edges, function(p) names(p$nodes))

  # Check for overlap
  all_genes <- unlist(pathway_genes)
  unique_genes <- unique(all_genes)

  # There should be shared genes (length of all_genes > unique_genes)
  expect_true(length(all_genes) > length(unique_genes))
})

test_that("multi-pathway networks have shared genes across pathways", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  # Check that genes appear in multiple pathways
  gene_id <- "gnaD"  # Should be in multiple ED pathways

  pathway_count <- 0
  for (pathway_id in names(potato@edges)) {
    pathway <- potato@edges[[pathway_id]]
    if (gene_id %in% names(pathway$nodes)) {
      pathway_count <- pathway_count + 1
    }
  }

  # Gene should appear in multiple pathways
  expect_true(pathway_count > 1)
})

test_that("compute_potato_hash works with multi-pathway networks", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  hash <- compute_potato_hash(potato)

  expect_type(hash, "character")
  expect_equal(nchar(hash), 32)  # MD5 hash length
})

test_that("multi-pathway scoring creates separate rows per pathway", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)

  # Copy ED network
  dir.create(file.path(temp_dir, "potatoes"), recursive = TRUE)
  file.copy(
    system.file("potatoes", "entner_doudoroff_network.json", package = "potato"),
    file.path(temp_dir, "potatoes", "entner_doudoroff_network.json")
  )

  sack <- create_sack(temp_dir)

  # Mock annotation results
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

  # Should have one row per pathway
  ed_scores <- sack@scores[sack@scores$potato == "entner_doudoroff_network", ]
  num_pathways <- length(potato@edges)

  expect_equal(nrow(ed_scores), num_pathways)
  expect_true(all(c("pathway", "pathway_name") %in% names(ed_scores)))

  unlink(temp_dir, recursive = TRUE)
})

test_that("pathway filtering works for multi-pathway networks", {
  potato <- load_potato(system.file("potatoes", "entner_doudoroff_network.json", package = "potato"))

  # Get first pathway ID
  pathway_id <- names(potato@edges)[1]

  # This should work when prepare_potato_for_plotting is called with pathway parameter
  # Testing the internal function behavior
  prep <- prepare_potato_for_plotting(potato, NULL, NULL, FALSE, pathway_id)

  expect_equal(length(prep$potato@edges), 1)
  expect_equal(names(prep$potato@edges)[1], pathway_id)
})
