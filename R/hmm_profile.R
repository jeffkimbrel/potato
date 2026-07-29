#' Create filtered HMM profile file from potato detection terms
#'
#' Extracts only the HMM profiles needed by potatoes from the configured
#' HMM files and creates a filtered HMM profile file (may be concatenated).
#'
#' @param sack PotatoSack object
#' @param potato_names Vector of potato names (NULL = all)
#'
#' @returns List with hmm_profile path and modified sack
#' @export
create_hmm_profile <- function(sack, potato_names = NULL) {

  # Filter potatoes
  if (is.null(potato_names)) {
    potatoes <- sack@potatoes
  } else {
    potatoes <- sack@potatoes[potato_names]
  }

  # Extract all unique HMM profile NAMEs across potatoes
  all_profiles <- character()
  for (potato in potatoes) {
    for (node in potato@nodes) {
      if (!is.null(node$databases$hmm)) {
        all_profiles <- c(all_profiles, node$databases$hmm)
      }
    }
  }
  all_profiles <- unique(all_profiles)

  if (length(all_profiles) == 0) {
    cli::cli_abort("No hmm terms found in selected potatoes")
  }

  cli::cli_alert_info("Creating HMM profile for {length(all_profiles)} profile{?s}")

  # Get hmm files from config
  hmm_config <- sack@config$databases$hmm
  hmm_files <- hmm_config$files

  if (is.null(hmm_files) || length(hmm_files) == 0) {
    cli::cli_abort("No hmm files configured in potato_config.yaml")
  }

  # Read all HMM files and extract needed profiles
  needed_profiles <- list()
  found_profiles <- character()

  for (hmm_file in hmm_files) {
    if (!file.exists(hmm_file)) {
      cli::cli_warn("HMM file not found: {.file {hmm_file}}")
      next
    }

    # Read HMM file and parse profiles
    lines <- readLines(hmm_file)
    current_profile <- list()
    current_name <- NULL
    in_profile <- FALSE

    for (i in seq_along(lines)) {
      line <- lines[i]

      # Start of new profile
      if (grepl("^HMMER", line)) {
        # Save previous profile if it was needed
        if (!is.null(current_name) && length(current_profile) > 0) {
          if (current_name %in% all_profiles) {
            needed_profiles[[current_name]] <- current_profile
            found_profiles <- c(found_profiles, current_name)
          }
        }
        # Start new profile
        current_profile <- c(line)
        current_name <- NULL
        in_profile <- TRUE
      } else if (grepl("^NAME\\s+", line)) {
        # Extract profile NAME
        current_name <- sub("^NAME\\s+", "", line)
        current_name <- trimws(current_name)
        current_profile <- c(current_profile, line)
      } else if (grepl("^//", line)) {
        # End of profile
        current_profile <- c(current_profile, line)
        # Save if needed
        if (!is.null(current_name) && current_name %in% all_profiles) {
          needed_profiles[[current_name]] <- current_profile
          found_profiles <- c(found_profiles, current_name)
        }
        # Reset for next profile
        current_profile <- list()
        current_name <- NULL
        in_profile <- FALSE
      } else if (in_profile) {
        # Part of current profile
        current_profile <- c(current_profile, line)
      }
    }

    # Don't forget the last profile if file doesn't end with //
    if (!is.null(current_name) && length(current_profile) > 0) {
      if (current_name %in% all_profiles) {
        needed_profiles[[current_name]] <- current_profile
        found_profiles <- c(found_profiles, current_name)
      }
    }
  }

  # Check for missing profiles
  missing_profiles <- setdiff(all_profiles, found_profiles)
  if (length(missing_profiles) > 0) {
    cli::cli_warn(c(
      "Missing HMM profiles:",
      "x" = "{.val {missing_profiles}}"
    ))
  }

  if (length(needed_profiles) == 0) {
    cli::cli_abort("No HMM profiles found in hmm files")
  }

  # Get or create annotation session timestamp
  if (is.null(sack@metadata$annotation_session)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    sack@metadata$annotation_session <- timestamp
  } else {
    timestamp <- sack@metadata$annotation_session
  }

  # Create results/annotations/{timestamp} directory if needed
  hmm_dir <- file.path(sack@sack_root, "results", "annotations", timestamp)
  if (!dir.exists(hmm_dir)) {
    dir.create(hmm_dir, recursive = TRUE)
  }

  # Write filtered HMM profile file (concatenated)
  filtered_hmm <- file.path(hmm_dir, "filtered_profiles.hmm")
  hmm_lines <- character()
  for (profile_name in names(needed_profiles)) {
    hmm_lines <- c(hmm_lines, needed_profiles[[profile_name]])
  }
  writeLines(hmm_lines, filtered_hmm)

  cli::cli_alert_success("Created filtered HMM with {length(needed_profiles)} profile{?s}")

  list(hmm_profile = filtered_hmm, sack = sack)
}
