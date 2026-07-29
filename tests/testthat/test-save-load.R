test_that("save and load RDS roundtrip works", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)

  # Load sack
  sack <- load_potato_sack(temp_dir)

  # Save to RDS
  rds_path <- file.path(temp_dir, "test_sack.rds")
  save_potato_sack(sack, rds_path)

  expect_true(file.exists(rds_path))

  # Load from RDS
  loaded_sack <- load_potato_sack(rds_path)

  expect_true(S7::S7_inherits(loaded_sack, PotatoSack))
  expect_equal(loaded_sack@sack_id, sack@sack_id)
  expect_equal(length(loaded_sack@potatoes), length(sack@potatoes))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("inspect_saved_sack shows metadata", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)

  sack <- load_potato_sack(temp_dir)
  rds_path <- file.path(temp_dir, "test_sack.rds")
  save_potato_sack(sack, rds_path)

  # Inspect should not error
  metadata <- inspect_saved_sack(rds_path)

  expect_type(metadata, "list")
  expect_true("created" %in% names(metadata))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})
