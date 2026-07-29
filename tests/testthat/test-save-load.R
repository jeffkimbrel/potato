test_that("saveRDS and readRDS roundtrip works", {
  temp_dir <- tempfile()
  initialize_potato_sack(temp_dir, copy_potatoes = TRUE)

  # Create sack
  sack <- create_sack(temp_dir)

  # Save to RDS using standard R function
  rds_path <- file.path(temp_dir, "test_sack.rds")
  saveRDS(sack, rds_path)

  expect_true(file.exists(rds_path))

  # Load from RDS using standard R function
  loaded_sack <- readRDS(rds_path)

  expect_true(S7::S7_inherits(loaded_sack, PotatoSack))
  expect_equal(loaded_sack@sack_id, sack@sack_id)
  expect_equal(length(loaded_sack@potatoes), length(sack@potatoes))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})
