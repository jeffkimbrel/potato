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
      cli::cli_alert_warning("Overwriting existing hmm results")
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

  cli::cli_alert_info("Preparing HMM annotation...")

  # Create filtered HMM profile from potato detection terms
  hmm_result <- create_hmm_profile(sack, potato_names)
  hmm_profile <- hmm_result$hmm_profile
  sack <- hmm_result$sack

  # Get annotation session directory (created by create_hmm_profile)
  annotation_dir <- file.path(sack@sack_root, "results", "annotations", sack@metadata$annotation_session)

  # Filter potatoes
  if (is.null(potato_names)) {
    potatoes <- sack@potatoes
  } else {
    potatoes <- sack@potatoes[potato_names]
  }

  # Convert potato S7 objects to raw list
  potato_data <- lapply(potatoes, function(p) {
    list(
      id = p@id,
      name = p@name,
      nodes = p@nodes
    )
  })

  # Extract genome info
  genome_names <- sapply(sack@genomes, function(g) g@short_name)
  genome_paths <- sapply(sack@genomes, function(g) g@file_path)

  cli::cli_alert_info("Running HMM on {length(genome_paths)} genome{?s}...")

  # STEP 1: Run hmmsearch commands in parallel (just execute, return raw output + command)
  run_hmm_cmd <- function(genome_path, genome_name, hmm_profile, conda_env) {
    # hmmsearch: --tblout gives parseable table format, --noali skips alignments
    temp_out <- tempfile(fileext = ".tbl")
    cmd <- sprintf("hmmsearch --tblout %s --noali %s %s > /dev/null",
                   shQuote(temp_out), shQuote(hmm_profile), shQuote(genome_path))

    if (!is.null(conda_env)) {
      cmd <- sprintf("conda run -n %s %s", conda_env, cmd)
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
        result <- run_hmm_cmd(genome_paths[i], genome_names[i], hmm_profile, conda_env)
        p()
        result
      }, .options = furrr::furrr_options(seed = TRUE))
    })
    future::plan(future::sequential)
  } else {
    cli::cli_progress_bar("Running HMM", total = length(genome_paths))
    results <- purrr::map(seq_along(genome_paths), function(i) {
      cli::cli_progress_update()
      run_hmm_cmd(genome_paths[i], genome_names[i], hmm_profile, conda_env)
    })
    cli::cli_progress_done()
  }

  # Extract outputs and commands
  raw_outputs <- purrr::map(results, ~.x$output)
  commands <- purrr::map_chr(results, ~.x$command)

  # Write raw outputs to files
  cli::cli_alert_info("Saving raw outputs...")
  for (i in seq_along(genome_names)) {
    output_file <- file.path(annotation_dir, paste0(genome_names[i], ".hmm.txt"))
    writeLines(raw_outputs[[i]], output_file)
  }

  # Write command log
  log_file <- file.path(annotation_dir, "hmm.log")
  log_lines <- paste0(genome_names, "\t", commands)
  writeLines(c("genome\tcommand", log_lines), log_file)
  cli::cli_alert_success("Saved outputs to {.path {annotation_dir}}")

  # STEP 2: Parse outputs sequentially
  cli::cli_alert_info("Parsing HMM results...")

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
    hmm_hits_to_tibble(parsed, potato_data)
  })

  # Add to sack results
  sack@results$hmm <- hmm_results

  cli::cli_alert_success("HMM annotation complete")

  sack
}


#' Convert HMM hits to tibble (internal)
#' @noRd
hmm_hits_to_tibble <- function(parsed_hits, potato_data) {
  # Build map of HMM profile name -> potato nodes
  profile_to_nodes <- list()
  for (potato_id in names(potato_data)) {
    potato <- potato_data[[potato_id]]
    for (node in potato$nodes) {
      if (!is.null(node$databases$hmm)) {
        for (profile_name in node$databases$hmm) {
          if (is.null(profile_to_nodes[[profile_name]])) {
            profile_to_nodes[[profile_name]] <- list()
          }
          profile_to_nodes[[profile_name]][[length(profile_to_nodes[[profile_name]]) + 1]] <- list(
            potato = potato_id,
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

    for (hit in hits_list) {
      for (node_info in nodes) {
        rows[[length(rows) + 1]] <- list(
          potato = node_info$potato,
          node_id = node_info$node_id,
          step = node_info$step,
          target = hit$target,
          query = profile_name,
          evalue = hit$evalue,
          score = hit$score,
          bias = hit$bias
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(tibble::tibble(
      potato = character(),
      node_id = character(),
      step = integer(),
      target = character(),
      query = character(),
      evalue = numeric(),
      score = numeric(),
      bias = numeric()
    ))
  }

  dplyr::bind_rows(rows)
}
