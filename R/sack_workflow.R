#' Load a potato sack context
#'
#' Automatically discovers and loads config, potatoes, and genomes from a
#' potato sack directory. If no path provided, searches upward from current
#' working directory.
#'
#' @param path Path to potato sack. If NULL, uses find_potato_sack() from cwd.
#' @param genomes_dir Optional path to genomes directory. If provided, overrides
#'   the genomes path in the sack. Use this to point to genomes stored elsewhere.
#'
#' @returns PotatoSack S7 object with config, potatoes, genomes, and metadata
#'
#' @importFrom utils head packageVersion
#' @export
load_potato_sack <- function(path = NULL, genomes_dir = NULL) {

  # Find sack root
  if (is.null(path)) {
    sack_root <- find_potato_sack()
    if (is.null(sack_root)) {
      stop("Not inside a potato sack. Provide path or run from within a sack.", call. = FALSE)
    }
  } else {
    sack_root <- normalizePath(path, mustWork = TRUE)
    if (!file.exists(file.path(sack_root, "potato_config.yaml"))) {
      stop("Not a potato sack: ", sack_root, "\n",
           "  (missing potato_config.yaml)", call. = FALSE)
    }
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Loading potato sack: ", sack_root)
  } else {
    cli::cli_h2("Loading potato sack: {.path {sack_root}}")
  }

  # Load config
  config_path <- file.path(sack_root, "potato_config.yaml")
  config <- load_potato_config(config_path)
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("  Config loaded")
  } else {
    cli::cli_alert_success("Config loaded")
  }

  # Load potatoes
  potatoes_dir <- file.path(sack_root,
                           if (!is.null(config$paths$potatoes)) config$paths$potatoes else "potatoes")
  potatoes <- NULL
  if (dir.exists(potatoes_dir)) {
    potatoes <- load_potatoes(potatoes_dir)
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("  Loaded ", length(potatoes), " potato(s)")
    } else {
      cli::cli_alert_success("Loaded {length(potatoes)} potato{?s}")
    }
  } else {
    if (!requireNamespace("cli", quietly = TRUE)) {
      warning("Potatoes directory not found: ", potatoes_dir)
    } else {
      cli::cli_alert_warning("Potatoes directory not found: {.path {potatoes_dir}}")
    }
  }

  # Discover genomes
  if (is.null(genomes_dir)) {
    # Use default genomes path from sack
    genomes_dir <- file.path(sack_root,
                            if (!is.null(config$paths$genomes)) config$paths$genomes else "genomes")
  } else {
    # Use provided override path
    genomes_dir <- normalizePath(genomes_dir, mustWork = FALSE)
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("  Using genomes from: ", genomes_dir)
    } else {
      cli::cli_alert_info("Using genomes from: {.path {genomes_dir}}")
    }
  }

  genomes <- NULL
  if (dir.exists(genomes_dir)) {
    genome_files <- list.files(genomes_dir, pattern = "\\.(faa|fasta|fa|gbk|gb|gbff)$",
                               full.names = TRUE)
    if (length(genome_files) > 0) {
      genomes <- prepare_genomes(genome_dir = genomes_dir)
      if (!requireNamespace("cli", quietly = TRUE)) {
        message("  Discovered ", length(genomes), " genome(s)")
      } else {
        cli::cli_alert_success("Discovered {length(genomes)} genome{?s}")
      }
    } else {
      if (!requireNamespace("cli", quietly = TRUE)) {
        message("  No genomes found in ", genomes_dir)
      } else {
        cli::cli_alert_info("No genomes found in {.path {genomes_dir}}")
      }
    }
  } else {
    if (!requireNamespace("cli", quietly = TRUE)) {
      warning("Genomes directory not found: ", genomes_dir)
    } else {
      cli::cli_alert_warning("Genomes directory not found: {.path {genomes_dir}}")
    }
  }

  # Generate sack ID (hash of sack_root + timestamp)
  sack_id <- digest::digest(paste(sack_root, Sys.time()), algo = "md5")
  sack_id <- substr(sack_id, 1, 8)  # Short ID

  # Create metadata
  metadata <- list(
    created = Sys.time(),
    potato_version = as.character(packageVersion("potato")),
    databases = lapply(config$databases, function(db) {
      db_meta <- list(type = db$type)

      # Calculate database hash for change detection
      if (db$type == "kofam" && !is.null(db$profiles_dir)) {
        db_meta$profiles_hash <- calculate_dir_hash(db$profiles_dir, pattern = "\\.hmm$")
        db_meta$ko_list_hash <- calculate_file_hash(db$ko_list)
      } else if (db$type == "blast" && !is.null(db$path)) {
        db_meta$db_hash <- calculate_file_hash(db$path)
      } else if ((db$type == "hmm" || db$type == "pfam") && !is.null(db$path)) {
        db_meta$hmm_hash <- calculate_file_hash(db$path)
      }

      db_meta
    }),
    potatoes = lapply(potatoes, function(p) {
      list(
        id = p@id,
        name = p@name,
        json_path = p@json_path,
        hash = digest::digest(p, algo = "md5")
      )
    })
  )

  # Create PotatoSack object
  PotatoSack(
    sack_id = sack_id,
    sack_root = sack_root,
    config = config,
    potatoes = if (is.null(potatoes)) list() else potatoes,
    genomes = if (is.null(genomes)) list() else genomes,
    metadata = metadata
  )
}


#' Annotate all genomes in a sack
#'
#' Runs annotation across all tools configured in the sack, automatically
#' discovering genomes and potatoes.
#'
#' @param sack PotatoSack object from load_potato_sack(), or path to sack directory
#' @param tools Character vector of which tools to run. NULL = all configured tools
#' @param cleanup Logical. Remove intermediate files after parsing (default TRUE)
#'
#' @returns Modified PotatoSack object with results and provenance added
#' @export
annotate_sack <- function(sack, tools = NULL, cleanup = TRUE) {

  # Allow passing path directly
  if (is.character(sack)) {
    sack <- load_potato_sack(sack)
  }

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("First argument must be a PotatoSack object or path", call. = FALSE)
  }

  if (is.null(sack@genomes) || length(sack@genomes) == 0) {
    stop("No genomes found in sack", call. = FALSE)
  }

  if (is.null(sack@potatoes) || length(sack@potatoes) == 0) {
    stop("No potatoes found in sack", call. = FALSE)
  }

  # Validate potatoes
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("\nValidating potatoes...")
  } else {
    cli::cli_h2("Validating potatoes")
  }

  validation_errors <- character(0)
  validation_warnings <- character(0)

  for (potato in sack@potatoes) {
    val <- validate_potato(potato)
    if (!val$valid) {
      validation_errors <- c(validation_errors,
                            paste0("  [", potato@id, "] ", val$errors))
    }
    if (length(val$warnings) > 0) {
      validation_warnings <- c(validation_warnings,
                              paste0("  [", potato@id, "] ", val$warnings))
    }
  }

  if (length(validation_errors) > 0) {
    stop("Potato validation failed:\n", paste(validation_errors, collapse = "\n"), call. = FALSE)
  }

  if (length(validation_warnings) > 0) {
    for (warn in validation_warnings) {
      if (!requireNamespace("cli", quietly = TRUE)) {
        message(warn)
      } else {
        cli::cli_alert_warning(warn)
      }
    }
  }

  # Validate database resources
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("\nValidating database resources...")
  } else {
    cli::cli_h2("Validating database resources")
  }

  # Check KOfam profiles
  kofam_dbs <- names(sack$config$databases)[sapply(sack$config$databases, function(db) db$type == "kofam")]
  if (length(kofam_dbs) > 0) {
    for (db_name in kofam_dbs) {
      db_config <- sack$config$databases[[db_name]]
      profiles_dir <- db_config$profiles_dir

      # Get all KO terms from potatoes
      all_kos <- character(0)
      for (potato in sack@potatoes) {
        kos <- get_detection_terms(potato, database_name = db_name)
        if (length(kos) == 0) {
          kos <- get_detection_terms(potato, tool_type = "ko")
        }
        all_kos <- c(all_kos, kos)
      }
      all_kos <- unique(all_kos)

      if (length(all_kos) > 0) {
        missing_kos <- character(0)
        for (ko in all_kos) {
          profile_path <- file.path(profiles_dir, paste0(ko, ".hmm"))
          if (!file.exists(profile_path)) {
            missing_kos <- c(missing_kos, ko)
          }
        }

        if (length(missing_kos) > 0) {
          if (!requireNamespace("cli", quietly = TRUE)) {
            message("  [", db_name, "] Missing HMM profiles for ", length(missing_kos), " KO(s): ",
                    paste(head(missing_kos, 5), collapse = ", "),
                    if (length(missing_kos) > 5) "..." else "")
          } else {
            cli::cli_alert_warning("[{db_name}] Missing HMM profiles for {length(missing_kos)} KO(s): {paste(head(missing_kos, 5), collapse = ', ')}{if (length(missing_kos) > 5) '...' else ''}")
          }
        } else {
          if (!requireNamespace("cli", quietly = TRUE)) {
            message("  [", db_name, "] All ", length(all_kos), " KO profiles found")
          } else {
            cli::cli_alert_success("[{db_name}] All {length(all_kos)} KO profile{?s} found")
          }
        }
      }
    }
  }

  # Check BLAST databases
  blast_dbs <- names(sack$config$databases)[sapply(sack$config$databases, function(db) db$type == "blast")]
  if (length(blast_dbs) > 0) {
    for (db_name in blast_dbs) {
      db_config <- sack$config$databases[[db_name]]

      # Get all BLAST terms from potatoes
      all_blast_terms <- character(0)
      for (potato in sack@potatoes) {
        terms <- get_detection_terms(potato, database_name = db_name)
        if (length(terms) == 0) {
          terms <- get_detection_terms(potato, tool_type = "blast_terms")
        }
        all_blast_terms <- c(all_blast_terms, terms)
      }
      all_blast_terms <- unique(all_blast_terms)

      if (!file.exists(db_config$path)) {
        if (!requireNamespace("cli", quietly = TRUE)) {
          validation_errors <- c(validation_errors,
                                paste0("  [", db_name, "] BLAST database not found: ", db_config$path))
        } else {
          cli::cli_alert_danger("[{db_name}] BLAST database not found: {db_config$path}")
          validation_errors <- c(validation_errors, "BLAST database missing")
        }
      } else if (length(all_blast_terms) > 0) {
        # Check if terms exist in database (read fasta headers)
        fasta_headers <- tryCatch({
          lines <- readLines(db_config$path, n = 50000)
          header_lines <- lines[grepl("^>", lines)]
          gsub("^>([^ ]+).*", "\\1", header_lines)
        }, error = function(e) character(0))

        missing_terms <- setdiff(all_blast_terms, fasta_headers)

        if (length(missing_terms) > 0) {
          if (!requireNamespace("cli", quietly = TRUE)) {
            message("  [", db_name, "] Missing ", length(missing_terms), " BLAST term(s): ",
                    paste(head(missing_terms, 5), collapse = ", "),
                    if (length(missing_terms) > 5) "..." else "")
          } else {
            cli::cli_alert_warning("[{db_name}] Missing {length(missing_terms)} BLAST term{?s}: {paste(head(missing_terms, 5), collapse = ', ')}{if (length(missing_terms) > 5) '...' else ''}")
          }
        } else {
          if (!requireNamespace("cli", quietly = TRUE)) {
            message("  [", db_name, "] All ", length(all_blast_terms), " BLAST term(s) found")
          } else {
            cli::cli_alert_success("[{db_name}] All {length(all_blast_terms)} BLAST term{?s} found in database")
          }
        }
      } else {
        if (!requireNamespace("cli", quietly = TRUE)) {
          message("  [", db_name, "] Database found: ", basename(db_config$path))
        } else {
          cli::cli_alert_success("[{db_name}] Database found: {basename(db_config$path)}")
        }
      }
    }
  }

  if (length(validation_errors) > 0) {
    stop("Database validation failed", call. = FALSE)
  }

  # Initialize table
  anno_table <- initialize_annotation_table(sack$genomes, sack$config)

  # Get available tools from config
  available_tools <- unique(sapply(sack$config$databases, function(db) db$type))

  if (!is.null(tools)) {
    available_tools <- intersect(available_tools, tools)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("\nRunning annotation with tools: ", paste(available_tools, collapse = ", "))
  } else {
    cli::cli_h2("Running annotation")
    cli::cli_alert_info("Tools: {paste(available_tools, collapse = ', ')}")
  }

  # Run each tool
  for (tool in available_tools) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("\n=== ", toupper(tool), " ===")
    } else {
      cli::cli_h3(toupper(tool))
    }

    if (tool == "kofam") {
      anno_table$kofam <- annotate_genomes_kofam(sack$genomes, sack$potatoes, sack$config, cleanup = cleanup)
    } else if (tool == "blast") {
      anno_table$blast <- annotate_genomes_blast(sack$genomes, sack$potatoes, sack$config, cleanup = cleanup)
    } else if (tool == "hmm") {
      # TODO: implement annotate_genomes_hmm()
      if (!requireNamespace("cli", quietly = TRUE)) {
        warning("HMM annotation not yet implemented")
      } else {
        cli::cli_alert_warning("HMM annotation not yet implemented")
      }
    } else if (tool == "pfam") {
      # TODO: implement annotate_genomes_pfam()
      if (!requireNamespace("cli", quietly = TRUE)) {
        warning("PFAM annotation not yet implemented")
      } else {
        cli::cli_alert_warning("PFAM annotation not yet implemented")
      }
    }
  }

  # Map to potatoes
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("\nMapping results to potatoes...")
  } else {
    cli::cli_h2("Mapping results to potatoes")
  }
  anno_table <- map_annotation_table(anno_table, sack@potatoes)

  # Add results to sack
  sack@results <- anno_table
  sack@completed_stages <- c(sack@completed_stages, "annotation")

  # Mark annotation as completed
  if (!requireNamespace("cli", quietly = TRUE)) {
    message("\nAnnotation complete!")
  } else {
    cli::cli_alert_success("Annotation complete!")
  }

  sack
}
