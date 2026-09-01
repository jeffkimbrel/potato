#' Run HMM annotation on all genomes
#'
#' @param sack PotatoSack object
#' @param potato_names Vector of potato names (NULL = all)
#' @param conda_env Optional conda environment name (defaults to config setting)
#' @param workers Number of parallel workers (defaults to config setting, 1 = sequential)
#' @param overwrite If FALSE (default), error if hmm results already exist
#'
#' @returns Modified PotatoSack with hmm results in @results$hmm
#' @export

run_hmm <- function(sack, potato_names = NULL, conda_env = NULL, workers = NULL, overwrite = FALSE) {

  # Initialize message collection
  messages_list <- list()

  # Get conda_env from config if not provided
  if (is.null(conda_env)) {
    conda_env <- sack@config$annotation$conda_env
  }

  # Get workers from config if not provided
  if (is.null(workers)) {
    workers <- sack@config$annotation$workers
    if (is.null(workers)) workers <- 1
  }

  # Check if hmm results already exist
  if (!is.null(sack@results) && "hmm" %in% names(sack@results)) {
    if (!overwrite) {
      cli::cli_abort(c(
        "hmm annotation results already exist in sack",
        "i" = "Use {.code overwrite = TRUE} to replace existing results"
      ))
    } else {
      msg <- "Overwriting existing hmm results"
      cli::cli_alert_warning(msg)
      messages_list[[length(messages_list) + 1]] <- list(type = "warning", message = msg)
      sack@results$hmm <- NULL
    }
  }

  # Check genomes exist
  if (length(sack@genomes) == 0) {
    cli::cli_abort(c(
      "No genomes in sack",
      "i" = "Use {.fn add_genomes} first"
    ))
  }

  # Get hmm config
  hmm_config <- sack@config$databases$hmm

  if (is.null(hmm_config)) {
    cli::cli_abort("No hmm database configured in potato_config.yaml")
  }

  msg <- "Preparing HMM annotation..."
  cli::cli_alert_info(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "info", message = msg)

  # Create filtered HMM profile from potato detection terms
  hmm_result <- create_hmm_profile(sack, potato_names)
  hmm_profile <- hmm_result$hmm_profile
  tc_values <- hmm_result$tc_values
  all_profiles <- hmm_result$all_profiles
  potatoes_with_hmm <- hmm_result$potatoes_with_hmm
  profile_content <- hmm_result$profile_content
  sack <- hmm_result$sack

  # Get annotation session directory (created by create_hmm_profile)
  annotation_dir <- file.path(sack@sack_root, "results", "annotations", sack@metadata$annotation_session)

  # Filter potatoes
  if (is.null(potato_names)) {
    potatoes <- sack@potatoes
  } else {
    potatoes <- sack@potatoes[potato_names]
  }

  # Compute potato hashes for version tracking
  potato_hashes <- get_potato_hashes(potatoes)

  # Convert potato S7 objects to raw list
  potato_data <- lapply(potatoes, function(p) {
    list(
      id = p@id,
      name = p@name,
      genes = p@genes
    )
  })

  # Extract genome info
  genome_names <- sapply(sack@genomes, function(g) g@short_name)
  genome_paths <- sapply(sack@genomes, function(g) g@file_path)

  cli::cli_alert_info("Running HMM on {length(genome_paths)} genome{?s}...")

  # Find conda executable if needed
  conda_cmd <- "conda"
  if (!is.null(conda_env)) {
    conda_cmd <- find_conda()
    if (conda_cmd == "") {
      cli::cli_abort(c(
        "{.code conda} not found",
        "i" = "Make sure conda is installed",
        "i" = "Searched in PATH, CONDA_EXE, and common install locations"
      ))
    }
  }

  # Capture tool version for provenance
  version_cmd <- if (!is.null(conda_env)) {
    sprintf("%s run -n %s hmmsearch -h 2>&1 | head -2", conda_cmd, conda_env)
  } else {
    "hmmsearch -h 2>&1 | head -2"
  }
  tool_version <- tryCatch({
    version_output <- suppressWarnings(system(version_cmd, intern = TRUE))
    # hmmsearch version is in first 2 lines
    version_line <- paste(version_output[1:2], collapse = " ")
    if (is.na(version_line) || nchar(version_line) == 0) "unknown" else version_line
  }, error = function(e) {
    "unknown"
  })

  # STEP 1: Run hmmsearch commands in parallel (just execute, return raw output + command)
  run_hmm_cmd <- function(genome_path, genome_name, hmm_profile, conda_cmd, conda_env) {
    # hmmsearch: --tblout gives parseable table format, --noali skips alignments
    temp_out <- tempfile(fileext = ".tbl")
    cmd <- sprintf("hmmsearch --tblout %s --noali %s %s > /dev/null",
                   shQuote(temp_out), shQuote(hmm_profile), shQuote(genome_path))

    if (!is.null(conda_env)) {
      cmd <- sprintf("%s run -n %s %s", conda_cmd, conda_env, cmd)
    }

    system(cmd)
    output <- readLines(temp_out)
    unlink(temp_out)

    list(output = output, command = cmd)
  }

  # Run in parallel
  if (workers > 1) {
    # Use cli-style progress bar
    progressr::handlers(progressr::handler_cli(
      format = "{cli::pb_spin} Running HMM [{cli::pb_current}/{cli::pb_total}] | {cli::pb_bar} {cli::pb_percent}"
    ))
    progressr::with_progress({
      p <- progressr::progressor(along = genome_paths)
      results <- furrr::future_map(seq_along(genome_paths), function(i) {
        result <- run_hmm_cmd(genome_paths[i], genome_names[i], hmm_profile, conda_cmd, conda_env)
        p()
        result
      }, .options = furrr::furrr_options(seed = TRUE))
    })
    future::plan(future::sequential)
  } else {
    results <- list()
    cli::cli_progress_bar("Running HMM", total = length(genome_paths))
    for (i in seq_along(genome_paths)) {
      results[[i]] <- run_hmm_cmd(genome_paths[i], genome_names[i], hmm_profile, conda_cmd, conda_env)
      cli::cli_progress_update()
    }
    cli::cli_progress_done()
  }

  # Extract outputs and commands
  raw_outputs <- purrr::map(results, ~.x$output)
  commands <- purrr::map_chr(results, ~.x$command)

  # Write raw outputs to files
  msg <- "Saving raw outputs..."
  cli::cli_alert_info(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "info", message = msg)
  for (i in seq_along(genome_names)) {
    output_file <- file.path(annotation_dir, paste0(genome_names[i], ".hmm.txt"))
    writeLines(raw_outputs[[i]], output_file)
  }

  # Write command log
  log_file <- file.path(annotation_dir, "hmm.log")
  log_lines <- paste0(genome_names, "\t", commands)
  writeLines(c("genome\tcommand", log_lines), log_file)
  msg <- paste0("Saved outputs to ", annotation_dir)
  cli::cli_alert_success(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "success", message = msg)

  # STEP 2: Parse outputs sequentially
  msg <- "Parsing HMM results..."
  cli::cli_alert_info(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "info", message = msg)

  hmm_results <- purrr::map(raw_outputs, function(output_lines) {
    # Parse hmmsearch --tblout format
    # Columns: target_name accession query_name accession evalue score bias ...
    parsed <- list()

    for (line in output_lines) {
      if (length(line) > 0 && !grepl("^#", line)) {
        # Split on whitespace, but need to handle multiple spaces
        fields <- strsplit(trimws(line), "\\s+")[[1]]
        if (length(fields) >= 7) {
          hit <- list(
            target = fields[1],      # Query gene (from genome)
            query = fields[3],       # HMM profile name
            evalue = as.numeric(fields[5]),
            score = as.numeric(fields[6]),
            bias = as.numeric(fields[7])
          )

          # Group by HMM profile name
          if (is.null(parsed[[hit$query]])) {
            parsed[[hit$query]] <- list()
          }
          parsed[[hit$query]][[length(parsed[[hit$query]]) + 1]] <- hit
        }
      }
    }

    # Match to potato nodes
    hmm_hits_to_tibble(parsed, potato_data, potato_hashes, tc_values)
  })

  # Add to sack results
  sack@results$hmm <- hmm_results

  # Add provenance tracking
  # Build command template with placeholders
  hmmsearch_template <- "hmmsearch --tblout {tblout} --noali {hmm_profile} {genome_path} > /dev/null"
  if (!is.null(conda_env)) {
    hmmsearch_template <- sprintf("%s run -n %s %s", conda_cmd, conda_env, hmmsearch_template)
  }

  # Store messages in provenance
  messages_tbl <- tibble::tibble(
    type = sapply(messages_list, function(x) x$type),
    message = sapply(messages_list, function(x) x$message)
  )

  sack@provenance$hmm <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    tool_version = tool_version,
    conda_env = conda_env,
    workers = workers,
    potatoes_requested = names(potatoes),
    potatoes_with_genes = potatoes_with_hmm,
    n_genomes = length(genome_paths),
    n_profiles = length(all_profiles),
    messages = messages_tbl,
    commands = list(
      profile_content = profile_content,
      hmm_profile = hmm_profile,
      hmmsearch_template = hmmsearch_template,
      genome_paths = genome_paths
    )
  )

  msg <- "HMM annotation complete"
  cli::cli_alert_success(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "success", message = msg)

  # Update stored messages
  sack@provenance$hmm$messages <- tibble::tibble(
    type = sapply(messages_list, function(x) x$type),
    message = sapply(messages_list, function(x) x$message)
  )

  sack
}


#' Convert HMM hits to tibble (internal)
#' @noRd
hmm_hits_to_tibble <- function(parsed_hits, potato_data, potato_hashes, tc_values) {
  # Build map of HMM profile name -> potato nodes
  profile_to_nodes <- list()
  for (potato_id in names(potato_data)) {
    potato <- potato_data[[potato_id]]
    for (node in potato$genes) {
      if (!is.null(node$databases$hmm)) {
        for (profile_name in node$databases$hmm) {
          if (is.null(profile_to_nodes[[profile_name]])) {
            profile_to_nodes[[profile_name]] <- list()
          }
          profile_to_nodes[[profile_name]][[length(profile_to_nodes[[profile_name]]) + 1]] <- list(
            potato = potato_id,
            potato_hash = potato_hashes[[potato_id]],
            node_id = node$id,
            step = node$step
          )
        }
      }
    }
  }

  rows <- list()
  for (profile_name in names(parsed_hits)) {
    hits_list <- parsed_hits[[profile_name]]
    nodes <- profile_to_nodes[[profile_name]]
    if (is.null(nodes)) next

    # Get TC value for this profile (NA if not present)
    tc_threshold <- tc_values[[profile_name]]
    if (is.null(tc_threshold)) tc_threshold <- NA_real_

    for (hit in hits_list) {
      for (node_info in nodes) {
        rows[[length(rows) + 1]] <- list(
          potato = node_info$potato,
          potato_hash = node_info$potato_hash,
          node_id = node_info$node_id,
          step = node_info$step,
          target = hit$target,
          query = profile_name,
          evalue = hit$evalue,
          score = hit$score,
          bias = hit$bias,
          tc_threshold = tc_threshold
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(tibble::tibble(
      potato = character(),
      potato_hash = character(),
      node_id = character(),
      step = integer(),
      target = character(),
      query = character(),
      evalue = numeric(),
      score = numeric(),
      bias = numeric(),
      tc_threshold = numeric()
    ))
  }

  dplyr::bind_rows(rows)
}
