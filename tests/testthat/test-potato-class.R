test_that("load_potato reads valid potato JSON", {
  potato_path <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  skip_if(potato_path == "", "Test potato not found")

  potato <- load_potato(potato_path)

  expect_true(S7::S7_inherits(potato, Potato))
  expect_equal(potato@id, "test_glycolysis")
  expect_true(length(potato@nodes) > 0)
})

test_that("validate_potato accepts valid potatoes", {
  potato_path <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  skip_if(potato_path == "", "Test potato not found")

  potato <- load_potato(potato_path)
  validation <- validate_potato(potato)

  expect_true(validation$valid)
  expect_equal(length(validation$errors), 0)
})

test_that("validate_potato rejects invalid database types", {
  potato_path <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  skip_if(potato_path == "", "Test potato not found")

  potato <- load_potato(potato_path)

  # Manually inject invalid database type
  potato@nodes[[1]]$databases <- list(invalid_db = list("K00001"))

  validation <- validate_potato(potato)

  expect_false(validation$valid)
  expect_true(any(grepl("invalid database type", validation$errors)))
})

test_that("get_enzyme_nodes filters correctly", {
  potato_path <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  skip_if(potato_path == "", "Test potato not found")

  potato <- load_potato(potato_path)
  enzyme_nodes <- get_enzyme_nodes(potato)

  expect_true(length(enzyme_nodes) > 0)
  expect_true(all(sapply(enzyme_nodes, function(n) n$type == "enzyme")))
})

test_that("get_detection_terms extracts KO IDs", {
  potato_path <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  skip_if(potato_path == "", "Test potato not found")

  potato <- load_potato(potato_path)
  ko_terms <- get_detection_terms(potato, "kofam")

  expect_true(length(ko_terms) > 0)
  expect_true(all(grepl("^K[0-9]{5}$", ko_terms)))
})
