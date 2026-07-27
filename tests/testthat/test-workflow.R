test_that("load_potato_sack works with test sack", {
  # Create minimal test sack
  sack_path <- tempfile("test_sack_")
  dir.create(sack_path)
  dir.create(file.path(sack_path, "potatoes"))
  dir.create(file.path(sack_path, "genomes"))

  # Copy test potato
  test_potato_src <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  file.copy(test_potato_src, file.path(sack_path, "potatoes", "test_glycolysis.json"))

  # Create minimal config
  config_lines <- c(
    "project_name: test",
    "databases:",
    "  kofam113:",
    "    type: kofam",
    "    profiles_dir: /fake/path",
    "    ko_list: /fake/path/ko_list"
  )
  writeLines(config_lines, file.path(sack_path, "potato_config.yaml"))

  # Load sack
  sack <- load_potato_sack(sack_path)

  expect_true(S7::S7_inherits(sack, PotatoSack))
  expect_equal(basename(sack@sack_root), basename(sack_path))
  expect_type(sack@config, "list")
  expect_type(sack@potatoes, "list")
  expect_equal(length(sack@potatoes), 1)
  expect_true("test_glycolysis" %in% names(sack@potatoes))

  # Cleanup
  unlink(sack_path, recursive = TRUE)
})


test_that("load_potato_sack discovers genomes from custom path", {
  # Create minimal test sack
  sack_path <- tempfile("test_sack_")
  dir.create(sack_path)
  dir.create(file.path(sack_path, "potatoes"))

  # Create genomes elsewhere
  genomes_path <- tempfile("genomes_")
  dir.create(genomes_path)

  # Create dummy genome file
  writeLines(c(">gene1", "MKKK"), file.path(genomes_path, "test.faa"))

  # Copy test potato
  test_potato_src <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  file.copy(test_potato_src, file.path(sack_path, "potatoes", "test_glycolysis.json"))

  # Create minimal config
  config_lines <- c(
    "project_name: test",
    "databases:",
    "  kofam113:",
    "    type: kofam",
    "    profiles_dir: /fake/path",
    "    ko_list: /fake/path/ko_list"
  )
  writeLines(config_lines, file.path(sack_path, "potato_config.yaml"))

  # Skip if jakomics not available (needed for genome validation)
  skip_if_not(exists("jakomics"), "jakomics not loaded")

  # Load sack with custom genomes path
  sack <- load_potato_sack(sack_path, genomes_dir = genomes_path)

  expect_true(S7::S7_inherits(sack, PotatoSack))
  expect_type(sack@genomes, "list")
  expect_equal(length(sack@genomes), 1)

  # Cleanup
  unlink(sack_path, recursive = TRUE)
  unlink(genomes_path, recursive = TRUE)
})


test_that("load_potato_config validates new database schema", {
  # Create temp config with new schema
  config_path <- tempfile(fileext = ".yaml")

  config_lines <- c(
    "project_name: test",
    "databases:",
    "  kofam113:",
    "    type: kofam",
    "    profiles_dir: /fake/kofam",
    "    ko_list: /fake/kofam/ko_list",
    "  gator_blast:",
    "    type: blast",
    "    path: /fake/gator.faa"
  )
  writeLines(config_lines, config_path)

  config <- load_potato_config(config_path)

  expect_type(config, "list")
  expect_type(config$databases, "list")
  expect_equal(length(config$databases), 2)
  expect_equal(config$databases$kofam113$type, "kofam")
  expect_equal(config$databases$gator_blast$type, "blast")

  # Cleanup
  unlink(config_path)
})


test_that("get_detection_terms works with new database schema", {
  # Create test potato with new schema
  potato_data <- list(
    id = "test_new",
    name = "Test New Schema",
    type = "pathway",
    nodes = list(
      list(
        id = "gene1",
        type = "enzyme",
        name = "Gene 1",
        databases = list(
          kofam113 = c("K00001", "K00002"),
          gator_blast = c("gene1_ref")
        )
      )
    ),
    edges = list(),
    tags = c("test"),
    source = "",
    notes = ""
  )

  # Write and load
  temp_file <- tempfile(fileext = ".json")
  jsonlite::write_json(potato_data, temp_file, auto_unbox = TRUE)
  potato <- load_potato(temp_file)

  # Test new schema
  kos <- get_detection_terms(potato, database_name = "kofam113")
  expect_type(kos, "character")
  expect_equal(length(kos), 2)
  expect_true("K00001" %in% kos)

  blast_terms <- get_detection_terms(potato, database_name = "gator_blast")
  expect_type(blast_terms, "character")
  expect_equal(length(blast_terms), 1)
  expect_equal(blast_terms, "gene1_ref")

  # Cleanup
  unlink(temp_file)
})


test_that("save and load PotatoSack works", {
  # Create minimal test sack
  sack_path <- tempfile("test_sack_")
  dir.create(sack_path)
  dir.create(file.path(sack_path, "potatoes"))
  dir.create(file.path(sack_path, "genomes"))

  # Copy test potato
  test_potato_src <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  file.copy(test_potato_src, file.path(sack_path, "potatoes", "test_glycolysis.json"))

  # Create minimal config
  config_lines <- c(
    "project_name: test",
    "databases:",
    "  kofam113:",
    "    type: kofam",
    "    profiles_dir: /fake/path",
    "    ko_list: /fake/path/ko_list"
  )
  writeLines(config_lines, file.path(sack_path, "potato_config.yaml"))

  # Load sack
  sack <- load_potato_sack(sack_path)

  # Save
  save_path <- tempfile(fileext = ".rds")
  save_potato_sack(sack, save_path)

  expect_true(file.exists(save_path))

  # Load
  sack2 <- load_saved_sack(save_path)

  expect_true(S7::S7_inherits(sack2, PotatoSack))
  expect_equal(sack2@sack_id, sack@sack_id)
  expect_equal(length(sack2@potatoes), 1)

  # Inspect
  meta <- inspect_saved_sack(save_path)
  expect_type(meta, "list")
  expect_true("created" %in% names(meta))

  # Cleanup
  unlink(sack_path, recursive = TRUE)
  unlink(save_path)
})


test_that("validate_potato checks new database field", {
  # Create test potato with new schema
  potato_data <- list(
    id = "test_validation",
    name = "Test Validation",
    type = "pathway",
    nodes = list(
      list(
        id = "gene1",
        type = "enzyme",
        name = "Gene 1",
        databases = list(
          kofam113 = c("K00001")
        )
      )
    ),
    edges = list(),
    tags = c("test"),
    source = "",
    notes = ""
  )

  temp_file <- tempfile(fileext = ".json")
  jsonlite::write_json(potato_data, temp_file, auto_unbox = TRUE)
  potato <- load_potato(temp_file)

  validation <- validate_potato(potato)

  expect_true(validation$valid)
  expect_equal(length(validation$errors), 0)

  # Cleanup
  unlink(temp_file)
})
