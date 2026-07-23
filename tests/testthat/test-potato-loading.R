test_that("load_potato works with test_glycolysis", {
  test_path <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  skip_if(test_path == "", "Test potato not found")

  potato <- load_potato(test_path)

  expect_true(S7::S7_inherits(potato, Potato))
  expect_equal(potato@id, "test_glycolysis")
  expect_equal(potato@name, "Test Glycolysis (3 steps)")
  expect_equal(length(potato@nodes), 3)
  expect_equal(length(potato@edges), 2)
  expect_true("test" %in% potato@tags)
})

test_that("load_test_potato convenience function works", {
  potato <- load_test_potato()

  expect_true(S7::S7_inherits(potato, Potato))
  expect_equal(potato@id, "test_glycolysis")
})

test_that("get_enzyme_nodes returns only enzyme nodes", {
  potato <- load_test_potato()
  enzymes <- get_enzyme_nodes(potato)

  expect_type(enzymes, "list")
  expect_equal(length(enzymes), 3)

  # Check all are enzyme type
  for (enzyme in enzymes) {
    expect_equal(enzyme$type, "enzyme")
  }
})

test_that("get_detection_terms extracts KO and EC terms", {
  potato <- load_test_potato()

  ko_terms <- get_detection_terms(potato, "ko")
  expect_type(ko_terms, "character")
  expect_true(length(ko_terms) > 0)
  expect_true("K00134" %in% ko_terms)  # gapA

  ec_terms <- get_detection_terms(potato, "ec")
  expect_type(ec_terms, "character")
  expect_true("1.2.1.12" %in% ec_terms)  # gapA EC
})

test_that("validate_potato accepts valid potatoes", {
  potato <- load_test_potato()
  validation <- validate_potato(potato)

  expect_type(validation, "list")
  expect_true(validation$valid)
  expect_equal(length(validation$errors), 0)
  expect_true("warnings" %in% names(validation))
})

test_that("build_potato_graph creates igraph object", {
  potato <- load_test_potato()
  g <- build_potato_graph(potato)

  expect_s3_class(g, "igraph")
  expect_equal(igraph::vcount(g), 3)
  expect_equal(igraph::ecount(g), 2)
  expect_true(igraph::is_dag(g))
})

test_that("print method works for Potato", {
  skip_if(identical(Sys.getenv("_R_CHECK_PACKAGE_NAME_"), "potato"))
  potato <- load_test_potato()

  expect_output(print(potato), "Potato: test_glycolysis")
  expect_output(print(potato), "Nodes: 3")
})

test_that("summary method works for Potato", {
  skip_if(identical(Sys.getenv("_R_CHECK_PACKAGE_NAME_"), "potato"))
  potato <- load_test_potato()

  output <- capture.output(summary(potato))
  output_text <- paste(output, collapse = "\n")

  expect_true(grepl("Test Glycolysis", output_text))
  expect_true(grepl("Enzyme nodes:", output_text))
  expect_true(grepl("valid", output_text))
})

test_that("load_potatoes loads multiple potatoes from directory", {
  potatoes_dir <- system.file("potatoes", package = "potato")
  skip_if(potatoes_dir == "", "Potatoes directory not found")

  potatoes <- load_potatoes(potatoes_dir)

  expect_type(potatoes, "list")
  expect_true(length(potatoes) >= 2)
  expect_true("test_glycolysis" %in% names(potatoes))
  expect_true("nitrogen_fixation" %in% names(potatoes))
})

test_that("load_potatoes filters by tags", {
  potatoes_dir <- system.file("potatoes", package = "potato")
  skip_if(potatoes_dir == "", "Potatoes directory not found")

  potatoes <- load_potatoes(potatoes_dir, tags = "test")

  expect_type(potatoes, "list")
  expect_true("test_glycolysis" %in% names(potatoes))

  # Check all have the test tag
  for (p in potatoes) {
    expect_true("test" %in% p@tags)
  }
})
