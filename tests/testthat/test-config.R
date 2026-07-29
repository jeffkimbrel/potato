test_that("load_potato_config reads YAML", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)

  config_path <- file.path(temp_dir, "potato_config.yaml")
  config <- load_potato_config(config_path)

  expect_s3_class(config, "potato_config")
  expect_true("databases" %in% names(config))
  expect_true("paths" %in% names(config))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("load_potato_config validates database types", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)

  config_path <- file.path(temp_dir, "potato_config.yaml")
  config <- load_potato_config(config_path)

  # Check kofam database config
  expect_true("kofam" %in% names(config$databases))
  expect_equal(config$databases$kofam$type, "kofam")

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("load_potato_config auto-finds config", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = FALSE)

  old_wd <- getwd()
  setwd(temp_dir)

  # Should auto-find
  config <- load_potato_config()
  expect_s3_class(config, "potato_config")

  # Restore
  setwd(old_wd)
  unlink(temp_dir, recursive = TRUE)
})
