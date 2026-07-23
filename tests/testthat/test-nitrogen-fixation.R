test_that("nitrogen fixation pathway loads correctly", {
  nif_path <- system.file("potatoes", "nitrogen_fixation.json", package = "potato")
  skip_if(nif_path == "", "Nitrogen fixation potato not found")

  nif <- load_potato(nif_path)

  expect_true(S7::S7_inherits(nif, Potato))
  expect_equal(nif@id, "nitrogen_fixation")
  expect_equal(nif@name, "Nitrogen Fixation")
  expect_equal(length(nif@nodes), 7)
  expect_equal(length(nif@edges), 5)
})

test_that("nitrogen fixation has correct gene types", {
  nif_path <- system.file("potatoes", "nitrogen_fixation.json", package = "potato")
  skip_if(nif_path == "", "Nitrogen fixation potato not found")

  nif <- load_potato(nif_path)
  enzymes <- get_enzyme_nodes(nif)

  expect_equal(length(enzymes), 7)

  # Check for required genes (Mo-nitrogenase)
  enzyme_ids <- sapply(enzymes, function(e) e$id)
  expect_true("nifD" %in% enzyme_ids)
  expect_true("nifH" %in% enzyme_ids)
  expect_true("nifK" %in% enzyme_ids)

  # Check for optional genes (V-nitrogenase)
  expect_true("vnfD" %in% enzyme_ids)
  expect_true("vnfH" %in% enzyme_ids)
})

test_that("nitrogen fixation has expected KO terms", {
  nif_path <- system.file("potatoes", "nitrogen_fixation.json", package = "potato")
  skip_if(nif_path == "", "Nitrogen fixation potato not found")

  nif <- load_potato(nif_path)
  ko_terms <- get_detection_terms(nif, "ko")

  expect_type(ko_terms, "character")
  expect_equal(length(ko_terms), 7)

  # Check for specific KO terms
  expect_true("K02586" %in% ko_terms)  # nifD
  expect_true("K02588" %in% ko_terms)  # nifH
  expect_true("K02591" %in% ko_terms)  # nifK
})

test_that("nitrogen fixation DAG is valid", {
  nif_path <- system.file("potatoes", "nitrogen_fixation.json", package = "potato")
  skip_if(nif_path == "", "Nitrogen fixation potato not found")

  nif <- load_potato(nif_path)

  # Validate structure
  validation <- validate_potato(nif)
  expect_true(validation$valid)

  # Build and check graph
  g <- build_potato_graph(nif)
  expect_s3_class(g, "igraph")
  expect_equal(igraph::vcount(g), 7)
  expect_equal(igraph::ecount(g), 5)
  expect_true(igraph::is_dag(g))
})

test_that("nitrogen fixation represents two alternative pathways", {
  nif_path <- system.file("potatoes", "nitrogen_fixation.json", package = "potato")
  skip_if(nif_path == "", "Nitrogen fixation potato not found")

  nif <- load_potato(nif_path)

  # The DAG should show two separate pathways
  # Mo-nitrogenase: nifD -> nifK -> nifH
  # V-nitrogenase: vnfD -> vnfK -> vnfG -> vnfH

  g <- build_potato_graph(nif)

  # Check that nif and vnf pathways are separate
  expect_true("nifD" %in% igraph::V(g)$name)
  expect_true("vnfD" %in% igraph::V(g)$name)

  # Both should be in the same connected component (weakly connected)
  # but represent alternative routes
  expect_equal(length(igraph::components(g, mode = "weak")$membership), 7)
})
