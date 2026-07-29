#' Save a PotatoSack object to disk
#'
#' Saves a PotatoSack object as an RDS file with all metadata, results, and
#' provenance information preserved.
#'
#' @param sack PotatoSack object to save
#' @param path Path to save the RDS file (default: saves to sack directory as "<sack_id>_potato_sack.rds")
#' @param compress Compression method (default: "gzip")
#'
#' @return Invisibly returns the file path
#' @export
save_potato_sack <- function(sack, path = NULL, compress = "gzip") {

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("Input must be a PotatoSack object", call. = FALSE)
  }

  # Generate default path if not provided (save in sack directory)
  if (is.null(path)) {
    path <- file.path(sack@sack_root, paste0(sack@sack_id, "_potato_sack.rds"))
  }

  # Normalize to full path for display
  full_path <- normalizePath(path, mustWork = FALSE)

  # Update save timestamp in metadata
  sack@metadata$last_saved <- Sys.time()

  # Save as RDS
  saveRDS(sack, path, compress = compress)

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Saved PotatoSack to: ", full_path)
  } else {
    cli::cli_alert_success("Saved PotatoSack to: {.path {full_path}}")
  }

  invisible(path)
}


#' Load a PotatoSack object from disk
#'
#' Loads a previously saved PotatoSack object from an RDS file.
#'
#' @param path Path to the RDS file
#' @param validate Logical. Validate the loaded object (default TRUE)
#'
#' @return PotatoSack object
#' @export
load_saved_sack <- function(path, validate = TRUE) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  # Load RDS
  sack <- readRDS(path)

  # Validate
  if (validate) {
    if (!S7::S7_inherits(sack, PotatoSack)) {
      stop("Loaded object is not a PotatoSack", call. = FALSE)
    }
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Loaded PotatoSack from: ", path)
    message("  Sack root: ", sack@sack_root)
    message("  Genomes: ", length(sack@genomes), ", Potatoes: ", length(sack@potatoes))
  } else {
    cli::cli_alert_success("Loaded PotatoSack from {.path {basename(path)}}")
    cli::cli_text("  Sack root: {.path {sack@sack_root}}")
    cli::cli_text("  {length(sack@genomes)} genome{?s}, {length(sack@potatoes)} potato{?es}")
  }

  sack
}


#' Print metadata about a saved PotatoSack
#'
#' Shows information about a saved sack without fully loading it.
#'
#' @param path Path to the RDS file
#'
#' @return Invisibly returns the metadata
#' @export
inspect_saved_sack <- function(path) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  sack <- readRDS(path)

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("File is not a PotatoSack object", call. = FALSE)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    cat("<Saved PotatoSack>\n")
    cat("  File:", path, "\n")
    cat("  ID:", sack@sack_id, "\n")
    cat("  Created:", format(sack@metadata$created), "\n")
    if (!is.null(sack@metadata$last_saved)) {
      cat("  Last saved:", format(sack@metadata$last_saved), "\n")
    }
    cat("  Potatoes:", length(sack@potatoes), "\n")
    cat("  Genomes:", length(sack@genomes), "\n")
    if (length(sack@completed_stages) > 0) {
      cat("  Stages:", paste(sack@completed_stages, collapse = ", "), "\n")
    }
  } else {
    cli::cli_h2("Saved PotatoSack: {.path {basename(path)}}")
    cli::cli_dl(c(
      "ID" = sack@sack_id,
      "Created" = format(sack@metadata$created),
      "Last saved" = if (!is.null(sack@metadata$last_saved)) format(sack@metadata$last_saved) else "N/A",
      "Potatoes" = length(sack@potatoes),
      "Genomes" = length(sack@genomes),
      "Stages" = if (length(sack@completed_stages) > 0) paste(sack@completed_stages, collapse = ", ") else "none"
    ))

    if (!is.null(sack@results)) {
      cli::cli_alert_info("Has annotation results: {nrow(sack@results)} genome{?s}")
    }

    if (!is.null(sack@scores)) {
      cli::cli_alert_info("Has pathway scores: {nrow(sack@scores)} result{?s}")
    }
  }

  invisible(sack@metadata)
}


