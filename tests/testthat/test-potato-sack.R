test_that("initialize_potato_sack creates project structure", {
  skip_on_cran()

  test_dir <- tempdir()
  project_name <- paste0("test_sack_", as.integer(Sys.time()))

  # Initialize sack
  project_path <- initialize_potato_sack(test_dir, project_name, copy_potatoes = TRUE)

  expect_true(dir.exists(project_path))

  # Check expected files
  expect_true(file.exists(file.path(project_path, paste0(project_name, ".Rproj"))))
  expect_true(file.exists(file.path(project_path, "potato_config.yaml")))
  expect_true(file.exists(file.path(project_path, "README.md")))
  expect_true(file.exists(file.path(project_path, ".gitignore")))

  # Check expected directories
  expect_true(dir.exists(file.path(project_path, "potatoes")))
  expect_true(dir.exists(file.path(project_path, "genomes")))
  expect_true(dir.exists(file.path(project_path, "results")))

  # Check potatoes were copied
  potato_files <- list.files(file.path(project_path, "potatoes"), pattern = "\\.json$")
  expect_true(length(potato_files) >= 2)

  # Cleanup
  unlink(project_path, recursive = TRUE)
})

test_that("initialize_potato_sack with copy_potatoes=FALSE", {
  skip_on_cran()

  test_dir <- tempdir()
  project_name <- paste0("test_sack_no_copy_", as.integer(Sys.time()))

  # Initialize without copying potatoes
  project_path <- initialize_potato_sack(test_dir, project_name, copy_potatoes = FALSE)

  expect_true(dir.exists(project_path))

  # Potatoes directory should exist but be empty (or nearly empty)
  potato_files <- list.files(file.path(project_path, "potatoes"), pattern = "\\.json$")
  expect_equal(length(potato_files), 0)

  # Cleanup
  unlink(project_path, recursive = TRUE)
})

test_that("initialize_potato_sack fails on existing directory", {
  skip_on_cran()

  test_dir <- tempdir()
  project_name <- paste0("test_sack_exists_", as.integer(Sys.time()))

  # Create directory first
  project_path <- file.path(test_dir, project_name)
  dir.create(project_path)

  # Should fail
  expect_error(
    initialize_potato_sack(test_dir, project_name),
    "already exists"
  )

  # Cleanup
  unlink(project_path, recursive = TRUE)
})

test_that("find_potato_sack finds project root", {
  skip_on_cran()

  test_dir <- tempdir()
  project_name <- paste0("test_sack_find_", as.integer(Sys.time()))

  # Initialize sack
  project_path <- initialize_potato_sack(test_dir, project_name, copy_potatoes = FALSE)

  # Should find it
  found <- find_potato_sack(project_path)
  expect_equal(normalizePath(found), normalizePath(project_path))

  # Should find from subdirectory
  subdir <- file.path(project_path, "genomes")
  found_from_subdir <- find_potato_sack(subdir)
  expect_equal(normalizePath(found_from_subdir), normalizePath(project_path))

  # Cleanup
  unlink(project_path, recursive = TRUE)
})

test_that("find_potato_sack returns NULL outside project", {
  # Should not find a potato_config.yaml here
  result <- find_potato_sack(tempdir())
  expect_null(result)
})

test_that("potato_config.yaml has expected structure", {
  skip_on_cran()

  test_dir <- tempdir()
  project_name <- paste0("test_sack_config_", as.integer(Sys.time()))

  # Initialize sack
  project_path <- initialize_potato_sack(test_dir, project_name, copy_potatoes = FALSE)

  config_path <- file.path(project_path, "potato_config.yaml")
  expect_true(file.exists(config_path))

  # Read config
  config_lines <- readLines(config_path)

  # Check for expected sections
  expect_true(any(grepl("project_name:", config_lines)))
  expect_true(any(grepl("tools:", config_lines)))
  expect_true(any(grepl("kofam:", config_lines)))
  expect_true(any(grepl("blast:", config_lines)))

  # Cleanup
  unlink(project_path, recursive = TRUE)
})
