test_that("initialize_potato_sack creates project structure", {
  temp_dir <- tempfile()

  # Create sack
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)

  # Check structure
  expect_true(dir.exists(temp_dir))
  expect_true(file.exists(file.path(temp_dir, "potato_config.yaml")))
  expect_true(dir.exists(file.path(temp_dir, "potatoes")))
  expect_true(dir.exists(file.path(temp_dir, "genomes")))
  expect_true(dir.exists(file.path(temp_dir, "results")))
  expect_true(dir.exists(file.path(temp_dir, ".claude")))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("create_sack constructs from directory", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)

  # Create sack
  sack <- create_sack(temp_dir)

  expect_true(S7::S7_inherits(sack, PotatoSack))
  expect_equal(sack@sack_id, basename(temp_dir))
  expect_true(length(sack@potatoes) > 0)  # Should have copied example potatoes
  expect_equal(length(sack@genomes), 0)  # No genomes yet

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("find_potato_sack finds config file", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)

  # Change to sack directory
  old_wd <- getwd()
  setwd(temp_dir)

  # Should find config
  found <- find_potato_sack()
  expect_equal(normalizePath(found), normalizePath(temp_dir))

  # Restore
  setwd(old_wd)
  unlink(temp_dir, recursive = TRUE)
})

test_that("find_potato_sack returns NULL outside sack", {
  # Create temp dir without config
  temp_dir <- tempfile()
  dir.create(temp_dir)

  old_wd <- getwd()
  setwd(temp_dir)

  # Should not find config
  found <- find_potato_sack()
  expect_null(found)

  # Restore
  setwd(old_wd)
  unlink(temp_dir, recursive = TRUE)
})
