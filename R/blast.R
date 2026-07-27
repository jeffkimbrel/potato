#' Run BLAST annotation on a genome
#'
#' Executes blastp via jakomics to search for genes against a BLAST database.
#'
#' @param genome_file jakomics FILE object
#' @param blast_terms Character vector of BLAST gene IDs to search for
#' @param config Potato config object
#' @param db_name Name of BLAST database in config (e.g., "gator_blast")
#' @param db_path Optional direct path to database (overrides db_name lookup)
#' @param output_dir Directory for output files (default: "results/blast/")
#'
#' @returns Data frame of BLAST hits
#' @keywords internal
run_blast_annotation <- function(genome_file, blast_terms, config,
                                 db_name = NULL, db_path = NULL,
                                 output_dir = "results/blast/") {

  # Check jakomics is available
  if (!exists("jakomics")) {
    stop("jakomics not loaded. Check conda environment 'potato' is activated.", call. = FALSE)
  }

  if (length(blast_terms) == 0) {
    message("No BLAST terms provided")
    return(data.frame())
  }

  # Determine database path
  if (!is.null(db_path)) {
    # Direct path override
    blast_db_path <- db_path
    blast_config <- list()  # No config for direct path
  } else if (!is.null(db_name)) {
    # Look up in config
    if (is.null(config$databases[[db_name]])) {
      stop("Database '", db_name, "' not found in config", call. = FALSE)
    }
    blast_config <- config$databases[[db_name]]
    if (blast_config$type != "blast") {
      stop("Database '", db_name, "' is not a BLAST database (type: ", blast_config$type, ")", call. = FALSE)
    }
    blast_db_path <- blast_config$path
  } else {
    stop("Must provide either db_name or db_path", call. = FALSE)
  }

  if (!file.exists(blast_db_path)) {
    stop("BLAST database not found: ", blast_db_path, call. = FALSE)
  }

  # Check blastp executable
  executable <- if (!is.null(blast_config$executable)) blast_config$executable else "blastp"
  if (!jakomics$utilities$check_executable(executable)) {
    stop("blastp executable not found: ", executable, ". Check conda environment is activated.", call. = FALSE)
  }

  # Get thresholds
  evalue_thresh <- if (!is.null(blast_config$thresholds$evalue)) {
    blast_config$thresholds$evalue
  } else {
    1e-7
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Running BLAST: ", genome_file$short_name, " vs ", basename(blast_db_path))
  } else {
    cli::cli_alert_info("Running BLAST: {genome_file$short_name} vs {.file {basename(blast_db_path)}}")
  }

  # Get conda environment path for PATH
  py_path <- reticulate::py_config()$python
  conda_env_bin <- dirname(py_path)

  # Run BLAST via jakomics
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(conda_env_bin, old_path, sep = ":"))

  tryCatch({
    results <- jakomics$blast$run_blast(
      type = "prot",
      q = genome_file$file_path,
      db = blast_db_path,
      e = evalue_thresh,
      make = FALSE,
      return_query_results = FALSE,  # Return by subject (database gene IDs)
      echo = FALSE
    )

    # Restore PATH
    Sys.setenv(PATH = old_path)

    # Filter to only blast_terms we're interested in
    filtered_results <- list()
    for (term in blast_terms) {
      if (term %in% names(results)) {
        filtered_results[[term]] <- results[[term]]
      }
    }

    if (length(filtered_results) == 0) {
      message("  No hits to target BLAST terms")
      return(data.frame())
    }

    # Convert to data frame
    rows_list <- list()
    for (subject_id in names(filtered_results)) {
      hits <- filtered_results[[subject_id]]

      for (hit in hits) {
        # Get bitscore threshold if specified
        bitscore_thresh <- if (!is.null(blast_config$thresholds$bitscore)) {
          blast_config$thresholds$bitscore
        } else {
          0
        }

        pident_thresh <- if (!is.null(blast_config$thresholds$pident)) {
          blast_config$thresholds$pident
        } else {
          0
        }

        passed <- (hit$bit_score >= bitscore_thresh &&
                   hit$percent >= pident_thresh &&
                   hit$eval <= evalue_thresh)

        rows_list <- c(rows_list, list(data.frame(
          genome = genome_file$short_name,
          query = hit$query,
          subject = hit$subject,
          pident = hit$percent,
          length = hit$alignment_length,
          evalue = hit$eval,
          bitscore = hit$bit_score,
          passed = passed,
          stringsAsFactors = FALSE
        )))
      }
    }

    result_df <- do.call(rbind, rows_list)
    message("  Found ", nrow(result_df), " BLAST hits (",
            sum(result_df$passed), " passed thresholds)")
    result_df

  }, error = function(e) {
    Sys.setenv(PATH = old_path)
    warning("BLAST annotation failed for ", genome_file$short_name, ": ", e$message)
    return(NULL)
  })
}


#' Annotate multiple genomes with BLAST
#'
#' @param genomes List of jakomics FILE objects
#' @param potatoes List of Potato objects
#' @param config Potato config object
#' @param cleanup Logical. Remove intermediate files after parsing (default TRUE)
#'
#' @returns List of results, one per genome
#' @export
annotate_genomes_blast <- function(genomes, potatoes, config, cleanup = TRUE) {

  # Extract all BLAST terms from potatoes
  if (inherits(potatoes, "Potato")) {
    potatoes <- list(potatoes)
  }

  # Find all blast databases in config
  blast_dbs <- names(config$databases)[sapply(config$databases, function(db) db$type == "blast")]

  if (length(blast_dbs) == 0) {
    stop("No BLAST databases configured in config", call. = FALSE)
  }

  # Collect BLAST terms for each database
  db_terms <- list()  # db_name -> list of terms

  for (potato in potatoes) {
    for (node in potato@nodes) {
      if (is.null(node$type) || node$type != "enzyme") next

      # New schema: check databases field
      if (!is.null(node$databases)) {
        for (db_name in names(node$databases)) {
          # Check if this is a blast database
          if (db_name %in% blast_dbs) {
            if (!db_name %in% names(db_terms)) {
              db_terms[[db_name]] <- list(
                db_name = db_name,
                db_path = NULL,
                terms = character(0)
              )
            }
            db_terms[[db_name]]$terms <- c(db_terms[[db_name]]$terms, node$databases[[db_name]])
          }
        }
      }

      # Legacy schema: check blast_db and blast_terms
      if (!is.null(node$blast_terms)) {
        db_name <- NULL
        db_path <- NULL

        if (!is.null(node$blast_db_path)) {
          db_path <- node$blast_db_path
          db_key <- db_path
        } else if (!is.null(node$blast_db)) {
          db_name <- node$blast_db
          db_key <- db_name
        } else {
          next
        }

        if (!db_key %in% names(db_terms)) {
          db_terms[[db_key]] <- list(
            db_name = db_name,
            db_path = db_path,
            terms = character(0)
          )
        }

        db_terms[[db_key]]$terms <- c(db_terms[[db_key]]$terms, node$blast_terms)
      }
    }
  }

  if (length(db_terms) == 0) {
    stop("No BLAST terms found in provided potatoes", call. = FALSE)
  }

  # Remove duplicates
  for (db_key in names(db_terms)) {
    db_terms[[db_key]]$terms <- unique(db_terms[[db_key]]$terms)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Found BLAST terms across ", length(db_terms), " database(s)")
  } else {
    cli::cli_alert_info("Found BLAST terms across {length(db_terms)} database{?s}")
  }

  # For now, only support one database
  if (length(db_terms) > 1) {
    stop("Multiple BLAST databases not yet supported in parallel annotation. ",
         "Run separately for each database.", call. = FALSE)
  }

  db_info <- db_terms[[1]]
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("  Database: ", db_info$db_name %||% basename(db_info$db_path),
            " (", length(db_info$terms), " terms)")
  } else {
    cli::cli_alert_info("Database: {db_info$db_name %||% basename(db_info$db_path)} ({length(db_info$terms)} term{?s})")
  }

  # Run annotation
  results <- lapply(genomes, function(genome) {
    run_blast_annotation(
      genome,
      blast_terms = db_info$terms,
      config = config,
      db_name = db_info$db_name,
      db_path = db_info$db_path
    )
  })

  names(results) <- sapply(genomes, function(g) g$short_name)

  # Note: BLAST via jakomics returns results directly, no intermediate files to clean up

  results
}
