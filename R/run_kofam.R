#' Run kofam annotation on all genomes
#'
#' @param sack PotatoSack object
#' @param potato_names Vector of potato names (NULL = all)
#' @param conda_env Optional conda environment name (defaults to config setting)
#' @param workers Number of parallel workers (defaults to config setting, 1 = sequential)
#' @param overwrite If FALSE (default), error if kofam results already exist
#'
#' @returns Modified PotatoSack with kofam results in @results$kofam
#' @export

run_kofam <- function(sack, potato_names = NULL, conda_env = NULL, workers = NULL, overwrite = FALSE) {

  # Get conda_env from config if not provided
  if (is.null(conda_env)) {
    conda_env <- sack@config$annotation$conda_env
  }

  # Get workers from config if not provided
  if (is.null(workers)) {
    workers <- sack@config$annotation$workers
    if (is.null(workers)) workers <- 1
  }

  # Check if kofam results already exist
  if (!is.null(sack@results) && "kofam" %in% names(sack@results)) {
    if (!overwrite) {
      cli::cli_abort(c(
        "kofam annotation results already exist in sack",
        "i" = "Use {.code overwrite = TRUE} to replace existing results"
      ))
    } else {
      cli::cli_alert_warning("Overwriting existing kofam results")
      sack@results$kofam <- NULL
    }
  }

  # Check genomes exist
  if (length(sack@genomes) == 0) {
    cli::cli_abort(c(
      "No genomes in sack",
      "i" = "Use {.fn add_genomes} first"
    ))
  }

  # Create .hal file for selected potatoes
  cli::cli_alert_info("Preparing kofam annotation...")
  hal_result <- create_kofam_hal(sack, potato_names)
  hal_path <- hal_result$hal_path
  hal_content <- hal_result$hal_content
  all_kos <- hal_result$all_kos
  potatoes_with_kofam <- hal_result$potatoes_with_kofam
  sack <- hal_result$sack

  # Get kofam config
  kofam_config <- sack@config$databases$kofam

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

  # Get annotation session directory (created by create_kofam_hal)
  annotation_dir <- file.path(sack@sack_root, "results", "annotations", sack@metadata$annotation_session)

  # Set up parallel if requested
  if (workers > 1 && requireNamespace("furrr", quietly = TRUE)) {
    cli::cli_alert_info("Setting up parallel processing with {workers} worker{?s}")
    future::plan(future::multisession, workers = workers)
  } else {
    workers <- 1
  }

  # Extract genome info
  genome_names <- sapply(sack@genomes, function(g) g@short_name)
  genome_paths <- sapply(sack@genomes, function(g) g@file_path)

  cli::cli_alert_info("Running kofam on {length(genome_paths)} genome{?s}...")

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
    sprintf("%s run -n %s exec_annotation --version 2>&1", conda_cmd, conda_env)
  } else {
    "exec_annotation --version 2>&1"
  }
  tool_version <- tryCatch({
    version_output <- suppressWarnings(system(version_cmd, intern = TRUE))
    # Take first non-empty line
    version_line <- version_output[nchar(version_output) > 0][1]
    if (is.na(version_line)) "unknown" else version_line
  }, error = function(e) {
    "unknown"
  })

  # STEP 1: Run kofam commands in parallel (just execute, return raw output + command)
  run_kofam_cmd <- function(genome_path, genome_name, hal_path, ko_list, conda_cmd, conda_env) {
    temp_dir <- file.path(tempdir(), paste0("kofam_", basename(genome_path)))
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

    cmd <- sprintf("exec_annotation --no-report-unannotated -k %s --tmp-dir %s '%s' --cpu 1 --profile %s -f detail-tsv",
                   ko_list, temp_dir, genome_path, hal_path)

    if (!is.null(conda_env)) {
      cmd <- sprintf("%s run -n %s %s", conda_cmd, conda_env, cmd)
    }

    output <- system(cmd, intern = TRUE)

    list(output = output, command = cmd)
  }

  # Run in parallel
  if (workers > 1) {
    # Use cli-style progress bar
    progressr::handlers(progressr::handler_cli(
      format = "{cli::pb_spin} Running kofam [{cli::pb_current}/{cli::pb_total}] | {cli::pb_bar} {cli::pb_percent}"
    ))
    progressr::with_progress({
      p <- progressr::progressor(along = genome_paths)
      results <- furrr::future_map(seq_along(genome_paths), function(i) {
        result <- run_kofam_cmd(genome_paths[i], genome_names[i], hal_path, kofam_config$ko_list, conda_cmd, conda_env)
        p()
        result
      }, .options = furrr::furrr_options(seed = TRUE))
    })
    future::plan(future::sequential)
  } else {
    cli::cli_progress_bar("Running kofam", total = length(genome_paths))
    results <- purrr::map(seq_along(genome_paths), function(i) {
      cli::cli_progress_update()
      run_kofam_cmd(genome_paths[i], genome_names[i], hal_path, kofam_config$ko_list, conda_cmd, conda_env)
    })
    cli::cli_progress_done()
  }

  # Extract outputs and commands
  raw_outputs <- purrr::map(results, ~.x$output)
  commands <- purrr::map_chr(results, ~.x$command)

  # Write raw outputs to files
  cli::cli_alert_info("Saving raw outputs...")
  for (i in seq_along(genome_names)) {
    output_file <- file.path(annotation_dir, paste0(genome_names[i], ".kofam.txt"))
    writeLines(raw_outputs[[i]], output_file)
  }

  # Write command log
  log_file <- file.path(annotation_dir, "kofam.log")
  log_lines <- paste0(genome_names, "\t", commands)
  writeLines(c("genome\tcommand", log_lines), log_file)
  cli::cli_alert_success("Saved outputs to {.path {annotation_dir}}")

  # STEP 2: Parse outputs sequentially using jakomics
  cli::cli_alert_info("Parsing kofam results...")
  kegg <- reticulate::import("jakomics.kegg", delay_load = FALSE)

  kofam_results <- purrr::map(raw_outputs, function(output_lines) {
    # Parse with jakomics - collect ALL hits regardless of pass/fail
    # Scoring phase will determine which hits to use
    hits <- list()
    for (line in output_lines) {
      if (length(line) > 0 && !grepl("^#", line)) {
        # Strip leading whitespace/asterisk (kofam's pass/fail indicator)
        clean_line <- sub("^[*\\s]+", "", line)
        if (nchar(clean_line) > 0) {
          hits[[length(hits) + 1]] <- kegg$KOFAM(clean_line, t_scale = 1.0, score_as_ratio = FALSE)
        }
      }
    }

    # Build dictionary of ALL hits (don't filter by passed status)
    # This replaces kegg$parse_kofam_hits which filters to passed=TRUE only
    parsed <- list()
    for (hit in hits) {
      ko <- hit$KO
      if (is.null(parsed[[ko]])) {
        parsed[[ko]] <- list()
      }
      parsed[[ko]][[length(parsed[[ko]]) + 1]] <- hit
    }

    # Match to potato nodes
    kofam_hits_to_tibble(parsed, potato_data, potato_hashes)
  })

  # Add to sack results
  sack@results$kofam <- kofam_results

  # Add provenance tracking
  # Build command template with placeholders
  cmd_template <- "exec_annotation --no-report-unannotated -k {ko_list} --tmp-dir {tmp_dir} '{genome_path}' --cpu 1 --profile {hal_path} -f detail-tsv"
  if (!is.null(conda_env)) {
    cmd_template <- sprintf("%s run -n %s %s", conda_cmd, conda_env, cmd_template)
  }

  sack@provenance$kofam <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    tool_version = tool_version,
    conda_env = conda_env,
    workers = workers,
    potatoes_requested = names(potatoes),
    potatoes_with_genes = potatoes_with_kofam,
    n_genomes = length(genome_paths),
    n_kos = length(all_kos),
    commands = list(
      hal_content = hal_content,
      hal_path = hal_path,
      ko_list = kofam_config$ko_list,
      exec_annotation_template = cmd_template,
      genome_paths = genome_paths
    )
  )

  cli::cli_alert_success("Kofam annotation complete")

  sack
}


#' Convert jakomics kofam hits to tibble (internal)
#' @noRd
kofam_hits_to_tibble <- function(parsed_hits, potato_data, potato_hashes) {
  # Build map of KO -> potato genes
  ko_to_genes <- list()
  for (potato_id in names(potato_data)) {
    potato <- potato_data[[potato_id]]
    for (gene in potato$genes) {
      if (!is.null(gene$databases$kofam)) {
        for (ko in gene$databases$kofam) {
          if (is.null(ko_to_genes[[ko]])) {
            ko_to_genes[[ko]] <- list()
          }
          ko_to_genes[[ko]][[length(ko_to_genes[[ko]]) + 1]] <- list(
            potato = potato_id,
            potato_hash = potato_hashes[[potato_id]],
            node_id = gene$id,
            step = gene$step
          )
        }
      }
    }
  }

  rows <- list()
  for (ko in names(parsed_hits)) {
    hits_list <- parsed_hits[[ko]]
    genes <- ko_to_genes[[ko]]
    if (is.null(genes)) next

    for (hit in hits_list) {
      for (node_info in genes) {
        rows[[length(rows) + 1]] <- list(
          potato = node_info$potato,
          potato_hash = node_info$potato_hash,
          node_id = node_info$node_id,
          step = node_info$step,
          gene_id = hit$gene,
          ko = ko,
          score = hit$score,
          evalue = hit$evalue,
          threshold = hit$threshold
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
      gene_id = character(),
      ko = character(),
      score = numeric(),
      evalue = numeric(),
      threshold = numeric()
    ))
  }

  dplyr::bind_rows(rows)
}
