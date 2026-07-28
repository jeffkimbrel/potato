#' Generate HAL file for KOfam annotation
#'
#' Creates a .hal file listing HMM profile paths for specified KO terms.
#' HAL files are named by the hash of their KO content, allowing reuse across
#' multiple genomes and potatoes with the same KO requirements.
#'
#' @param ko_terms Character vector of KO IDs (e.g., c("K00001", "K00002"))
#' @param config Potato config object from load_potato_config()
#' @param output_dir Directory to write .hal file (default: "results/kofam/")
#' @param overwrite Logical. If FALSE, reuses existing .hal file with same hash
#'
#' @returns Path to the generated .hal file
#'
#' @details
#' The .hal file format is:
#' ```
#' /path/to/profiles/K00001.hmm
#' /path/to/profiles/K00002.hmm
#' ```
#'
#' File is named: `{hash}.hal` where hash is MD5 of sorted KO terms.
#'
#' @importFrom utils head
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' config <- load_potato_config()
#' potato <- load_test_potato()
#' ko_terms <- get_detection_terms(potato, "ko")
#' hal_file <- generate_hal_file(ko_terms, config)
#' }
generate_hal_file <- function(ko_terms, config, output_dir = "results/kofam/",
                              overwrite = FALSE) {

  if (length(ko_terms) == 0) {
    stop("No KO terms provided", call. = FALSE)
  }

  # Find first kofam database in config
  kofam_dbs <- names(config$databases)[sapply(config$databases, function(db) db$type == "kofam")]
  if (length(kofam_dbs) == 0) {
    stop("No kofam databases configured", call. = FALSE)
  }

  kofam_config <- config$databases[[kofam_dbs[1]]]

  profiles_dir <- kofam_config$profiles_dir
  if (!dir.exists(profiles_dir)) {
    stop("kofam profiles_dir not found: ", profiles_dir, call. = FALSE)
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Sort KO terms for consistent hashing
  ko_terms <- sort(unique(ko_terms))

  # Generate hash for filename
  ko_string <- paste(ko_terms, collapse = ",")
  hash <- digest::digest(ko_string, algo = "md5")
  hal_path <- file.path(output_dir, paste0(hash, ".hal"))

  # Check if exists and skip if not overwriting
  if (file.exists(hal_path) && !overwrite) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("Using existing .hal file: ", basename(hal_path))
    } else {
      cli::cli_alert_info("Using existing .hal file: {.file {basename(hal_path)}}")
    }
    return(hal_path)
  }

  # Build profile paths
  profile_paths <- character(length(ko_terms))
  missing_kos <- character(0)

  for (i in seq_along(ko_terms)) {
    ko <- ko_terms[i]
    profile_path <- file.path(profiles_dir, paste0(ko, ".hmm"))

    if (!file.exists(profile_path)) {
      missing_kos <- c(missing_kos, ko)
    } else {
      profile_paths[i] <- profile_path
    }
  }

  # Warn about missing profiles
  if (length(missing_kos) > 0) {
    warning("Missing HMM profiles for ", length(missing_kos), " KO terms: ",
            paste(head(missing_kos, 10), collapse = ", "),
            if (length(missing_kos) > 10) " ..." else "")
  }

  # Remove missing paths
  profile_paths <- profile_paths[nzchar(profile_paths)]

  if (length(profile_paths) == 0) {
    stop("No valid HMM profiles found for provided KO terms", call. = FALSE)
  }

  # Write .hal file
  writeLines(profile_paths, hal_path)
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Generated .hal file: ", basename(hal_path), " (", length(profile_paths), " profiles)")
  } else {
    cli::cli_alert_success("Generated .hal file: {.file {basename(hal_path)}} ({length(profile_paths)} profiles)")
  }

  hal_path
}


#' Run KOfam annotation on a genome
#'
#' Executes exec_annotation via jakomics to annotate a genome against a .hal file
#' of KOfam HMM profiles.
#'
#' @param genome_file jakomics FILE object (from prepare_genomes())
#' @param hal_file Path to .hal file (from generate_hal_file())
#' @param config Potato config object
#' @param output_dir Directory for output files (default: "results/kofam/")
#' @param db_name Name of database being used (for tracking in results)
#'
#' @returns Data frame of KOfam hits with columns: gene_id, ko, score, evalue, etc.
#'
#' @keywords internal
run_kofam_annotation <- function(genome_file, hal_file, config,
                                 output_dir = "results/kofam/",
                                 db_name = NULL) {

  # Check jakomics is available
  if (!exists("jakomics")) {
    stop("jakomics not loaded. Check conda environment 'potato' is activated.", call. = FALSE)
  }

  if (!file.exists(hal_file)) {
    stop("HAL file not found: ", hal_file, call. = FALSE)
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Find first kofam database in config
  kofam_dbs <- names(config$databases)[sapply(config$databases, function(db) db$type == "kofam")]
  if (length(kofam_dbs) == 0) {
    stop("No kofam databases configured", call. = FALSE)
  }

  kofam_config <- config$databases[[kofam_dbs[1]]]

  # Check executable is available
  executable <- if (is.null(kofam_config$executable)) "exec_annotation" else kofam_config$executable
  if (!jakomics$utilities$check_executable(executable)) {
    stop("kofam executable not found: ", executable, "\n",
         "  Check that it's in your PATH or conda environment is activated", call. = FALSE)
  }

  # Set up paths
  ko_list <- kofam_config$ko_list
  score_ratio <- if (is.null(kofam_config$thresholds$score_ratio)) 1.0 else kofam_config$thresholds$score_ratio

  # if (!requireNamespace("cli", quietly = TRUE)) {
  #   message("Annotating: ", genome_file$short_name)
  # } else {
  #   cli::cli_alert_info("Annotating: {genome_file$short_name}")
  # }

  # Create output file path
  output_file <- file.path(output_dir, paste0(genome_file$short_name, "_kofam.tsv"))

  # Run kofam via system command with output redirect
  # exec_annotation output doesn't capture well through subprocess pipes
  # Use processx to set PATH from reticulate's Python environment

  # Get conda environment path from reticulate
  py_path <- reticulate::py_config()$python
  conda_env_bin <- dirname(py_path)

  # Build command
  command <- paste0(
    executable, ' --no-report-unannotated ',
    '-k ', ko_list, ' ',
    '--tmp-dir ', output_dir, ' ',
    '"', genome_file$file_path, '" ',
    '-T ', if (score_ratio == 1.0) 1 else score_ratio, ' ',
    '--cpu 1 ',
    '--profile ', hal_file, ' ',
    '-f detail-tsv ',
    '> ', output_file
  )

  # Run with conda env in PATH
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(conda_env_bin, old_path, sep = ":"))

  exit_code <- system(command, ignore.stdout = FALSE, ignore.stderr = FALSE)

  # Restore PATH
  Sys.setenv(PATH = old_path)

  if (exit_code != 0) {
    # Store as attribute for collection
    result <- data.frame()
    attr(result, "message") <- list(
      level = "error",
      genome = genome_file$short_name,
      message = paste("KOfam annotation failed (exit code", exit_code, ")")
    )
    return(result)
  }

  # Parse output file
  if (!file.exists(output_file)) {
    result <- data.frame()
    attr(result, "message") <- list(
      level = "error",
      genome = genome_file$short_name,
      message = paste("KOfam output file not created:", output_file)
    )
    return(result)
  }

  # Read and parse results
  tryCatch({
    lines <- readLines(output_file)

    # Filter out comment lines
    data_lines <- lines[!grepl("^\\s*#", lines) & nzchar(lines)]

    if (length(data_lines) == 0) {
      result <- data.frame()
      attr(result, "message") <- list(
        level = "info",
        genome = genome_file$short_name,
        message = "No KO hits found"
      )
      return(result)
    }

    # Parse into data frame
    results <- do.call(rbind, lapply(data_lines, function(line) {
      parts <- strsplit(line, "\t")[[1]]

      if (length(parts) < 7) return(NULL)

      threshold <- as.numeric(parts[4])
      score <- as.numeric(parts[5])
      passed <- grepl("^\\*", parts[1])  # * prefix means passed threshold

      # Generate threshold message if failed
      threshold_message <- NA_character_
      if (!passed && !is.na(threshold) && !is.na(score)) {
        threshold_message <- sprintf("score %.1f < threshold %.1f", score, threshold)
      }

      data.frame(
        genome = genome_file$short_name,
        database = if (!is.null(db_name)) db_name else NA_character_,
        gene_id = parts[2],
        ko = parts[3],
        threshold = threshold,
        score = score,
        evalue = as.numeric(parts[6]),
        description = gsub('^"|"$', '', parts[7]),
        passed = passed,
        threshold_message = threshold_message,
        stringsAsFactors = FALSE
      )
    }))

    # message("  Found ", nrow(results), " KO hits (", sum(results$passed), " passed threshold)")
    results

  }, error = function(e) {
    result <- data.frame()
    attr(result, "message") <- list(
      level = "error",
      genome = genome_file$short_name,
      message = paste("Failed to parse KOfam results:", e$message)
    )
    return(result)
  })
}


#' Annotate multiple genomes with KOfam
#'
#' Run KOfam annotation across multiple genomes sequentially.
#'
#' @param genomes List of jakomics FILE objects
#' @param potatoes List of Potato objects (or single Potato)
#' @param config Potato config object
#' @param cleanup Logical. Remove intermediate files (.hal, .tsv) after parsing (default TRUE)
#'
#' @returns List of results, one per genome
#'
#' @export
annotate_genomes_kofam <- function(genomes, potatoes, config, cleanup = TRUE) {

  # Extract all KO terms from potatoes
  if (inherits(potatoes, "Potato")) {
    potatoes <- list(potatoes)
  }

  # Find all kofam databases in config
  kofam_dbs <- names(config$databases)[sapply(config$databases, function(db) db$type == "kofam")]

  if (length(kofam_dbs) == 0) {
    stop("No kofam databases configured in config", call. = FALSE)
  }

  # For now, use the first kofam database
  # TODO: Support multiple kofam databases
  kofam_db <- kofam_dbs[1]
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Using kofam database: ", kofam_db)
  } else {
    cli::cli_alert_info("Using kofam database: {kofam_db}")
  }

  all_ko_terms <- character(0)
  for (potato in potatoes) {
    # Try new schema first
    ko_terms <- get_detection_terms(potato, database_name = kofam_db)
    # Fall back to legacy
    if (length(ko_terms) == 0) {
      ko_terms <- get_detection_terms(potato, tool_type = "ko")
    }
    all_ko_terms <- c(all_ko_terms, ko_terms)
  }

  all_ko_terms <- unique(all_ko_terms)

  if (length(all_ko_terms) == 0) {
    stop("No KO terms found in provided potatoes for database ", kofam_db, call. = FALSE)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Found ", length(all_ko_terms), " unique KO terms across ", length(potatoes), " potato(es)")
  } else {
    cli::cli_alert_info("Found {length(all_ko_terms)} unique KO term{?s} across {length(potatoes)} potato{?es}")
  }

  # Generate HAL file (shared across all genomes)
  hal_file <- generate_hal_file(all_ko_terms, config)

  # Run annotation on each genome with progress bar
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required. Install with: install.packages('purrr')", call. = FALSE)
  }

  output_files <- character(0)

  # Collect messages during annotation
  collected_messages <- list()

  results <- purrr::map(genomes, function(genome) {
    result <- run_kofam_annotation(genome, hal_file, config, db_name = kofam_db)

    # Collect message if present
    msg <- attr(result, "message")
    if (!is.null(msg)) {
      msg$stage <- "kofam"
      collected_messages <<- c(collected_messages, list(msg))
    }

    # Track output file for cleanup
    output_file <- file.path("results/kofam", paste0(genome$short_name, "_kofam.tsv"))
    if (file.exists(output_file)) {
      output_files <<- c(output_files, output_file)
    }
    result
  }, .progress = paste0("KOfam [", kofam_db, "]"))

  # Name results by genome
  names(results) <- sapply(genomes, function(g) g$short_name)

  # Attach messages to results for collection at workflow level
  attr(results, "messages") <- collected_messages

  # Cleanup intermediate files
  if (cleanup) {
    # Remove HAL file
    if (file.exists(hal_file)) {
      unlink(hal_file)
    }
    # Remove TSV output files
    for (f in output_files) {
      if (file.exists(f)) {
        unlink(f)
      }
    }
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("Cleaned up intermediate files")
    } else {
      cli::cli_alert_success("Cleaned up intermediate files")
    }
  }

  results
}
