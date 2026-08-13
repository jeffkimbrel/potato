#' Create .hal file for kofam from potato KO list
#'
#' @param sack PotatoSack object
#' @param potato_names Vector of potato names (NULL = all)
#'
#' @returns List with hal_path and modified sack
#' @export
create_kofam_hal <- function(sack, potato_names = NULL) {

  # Filter potatoes
  if (is.null(potato_names)) {
    potatoes <- sack@potatoes
  } else {
    # Validate potato names exist
    available_names <- names(sack@potatoes)
    invalid_names <- setdiff(potato_names, available_names)

    if (length(invalid_names) > 0) {
      cli::cli_abort(c(
        "Invalid potato names: {.val {invalid_names}}",
        "i" = "Available potatoes: {.val {available_names}}"
      ))
    }

    potatoes <- sack@potatoes[potato_names]

  }

  # Extract all unique KO IDs across potatoes
  all_kos <- character()
  for (potato in potatoes) {
    for (gene in potato@genes) {
      if (!is.null(gene$databases$kofam)) {
        all_kos <- c(all_kos, gene$databases$kofam)
      }
    }
  }
  all_kos <- unique(all_kos)

  if (length(all_kos) == 0) {
    cli::cli_abort("No kofam KOs found in selected potatoes")
  }

  cli::cli_alert_info("Creating .hal file for {length(all_kos)} KO{?s}")

  # Get kofam profiles directory from config
  profiles_dir <- sack@config$databases$kofam$profiles_dir

  # Build paths to .hmm files
  hmm_paths <- file.path(profiles_dir, paste0(all_kos, ".hmm"))

  # Check that HMM files exist
  missing <- hmm_paths[!file.exists(hmm_paths)]
  if (length(missing) > 0) {
    cli::cli_warn(c(
      "Missing HMM files:",
      "x" = "{.file {basename(missing)}}"
    ))
    hmm_paths <- hmm_paths[file.exists(hmm_paths)]
  }

  if (length(hmm_paths) == 0) {
    cli::cli_abort("No valid HMM files found")
  }

  # Create content for .hal file
  hal_content <- paste(hmm_paths, collapse = "\n")

  # Generate MD5 hash of content
  hal_hash <- digest::digest(hal_content, algo = "md5", serialize = FALSE)

  # Get or create annotation session timestamp
  if (is.null(sack@metadata$annotation_session)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    sack@metadata$annotation_session <- timestamp
  } else {
    timestamp <- sack@metadata$annotation_session
  }

  # Create results/annotations/{timestamp} directory if needed
  hal_dir <- file.path(sack@sack_root, "results", "annotations", timestamp)
  if (!dir.exists(hal_dir)) {
    dir.create(hal_dir, recursive = TRUE)
  }

  # Write .hal file
  hal_path <- file.path(hal_dir, paste0(hal_hash, ".hal"))

  if (!file.exists(hal_path)) {
    writeLines(hal_content, hal_path)
    cli::cli_alert_success("Created {.file {basename(hal_path)}}")
  } else {
    cli::cli_alert_info("Using existing {.file {basename(hal_path)}}")
  }

  list(hal_path = hal_path, sack = sack)
}
