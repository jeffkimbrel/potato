test_that("sack_messages returns empty for new sack", {
  # Create minimal sack structure for testing
  sack <- PotatoSack(
    sack_id = "test",
    sack_root = tempdir(),
    config = list(databases = list()),
    messages = list()
  )

  msgs <- sack_messages(sack, as_dataframe = TRUE)
  expect_s3_class(msgs, "data.frame")
  expect_equal(nrow(msgs), 0)
})

test_that("add_sack_message trims whitespace", {
  sack <- PotatoSack(
    sack_id = "test",
    sack_root = tempdir(),
    config = list(databases = list()),
    messages = list()
  )

  sack <- add_sack_message(sack, level = "info", stage = "test",
                          message = "  Test message with spaces  ")

  msgs <- sack_messages(sack, as_dataframe = TRUE)
  expect_equal(msgs$message[1], "Test message with spaces")
})

test_that("sack_messages filters by level", {
  sack <- PotatoSack(
    sack_id = "test",
    sack_root = tempdir(),
    config = list(databases = list()),
    messages = list()
  )

  sack <- add_sack_message(sack, level = "info", stage = "test", message = "Info message")
  sack <- add_sack_message(sack, level = "warning", stage = "test", message = "Warning message")
  sack <- add_sack_message(sack, level = "error", stage = "test", message = "Error message")

  all_msgs <- sack_messages(sack)
  expect_equal(nrow(all_msgs), 3)

  warnings_only <- sack_messages(sack, level = "warning")
  expect_equal(nrow(warnings_only), 1)
  expect_equal(warnings_only$level[1], "warning")
})

test_that("sack_messages filters by stage", {
  sack <- PotatoSack(
    sack_id = "test",
    sack_root = tempdir(),
    config = list(databases = list()),
    messages = list()
  )

  sack <- add_sack_message(sack, level = "info", stage = "annotation", message = "Annotation message")
  sack <- add_sack_message(sack, level = "info", stage = "scoring", message = "Scoring message")

  anno_msgs <- sack_messages(sack, stage = "annotation")
  expect_equal(nrow(anno_msgs), 1)
  expect_equal(anno_msgs$stage[1], "annotation")
})

test_that("sack_msg both prints and collects", {
  sack <- PotatoSack(
    sack_id = "test",
    sack_root = tempdir(),
    config = list(databases = list()),
    messages = list()
  )

  # With verbose = FALSE, should only collect
  sack <- sack_msg(sack, "Test header", level = "header", stage = "test", verbose = FALSE)

  msgs <- sack_messages(sack)
  expect_equal(nrow(msgs), 1)
  expect_equal(msgs$level[1], "info")  # header maps to info
  expect_equal(msgs$message[1], "Test header")
})
