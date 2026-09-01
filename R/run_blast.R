#' Run BLAST annotation on all genomes
#'
#' @param sack PotatoSack object
#' @param potato_names Vector of potato names (NULL = all)
#' @param conda_env Optional conda environment name (defaults to config setting)
#' @param workers Number of parallel workers (defaults to config setting, 1 = sequential)
#' @param overwrite If FALSE (default), error if blast results already exist
#'
#' @returns Modified PotatoSack with blast results in @results$blast
#' @export

run_blast <- function(sack, potato_names = NULL, conda_env = NULL, workers = NULL, overwrite = FALSE) {

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

  # Check if blast results already exist
  if (!is.null(sack@results) && "blast" %in% names(sack@results)) {
    if (!overwrite) {
      cli::cli_abort(c(
        "blast annotation results already exist in sack",
        "i" = "Use {.code overwrite = TRUE} to replace existing results"
      ))
    } else {
      msg <- "Overwriting existing blast results"
      cli::cli_alert_warning(msg)
      messages_list[[length(messages_list) + 1]] <- list(type = "warning", message = msg)
      sack@results$blast <- NULL
    }
  }

  # Check genomes exist
  if (length(sack@genomes) == 0) {
    cli::cli_abort(c(
      "No genomes in sack",
      "i" = "Use {.fn add_genomes} first"
    ))
  }

  # Get blast config
  blast_config <- sack@config$databases$blast

  if (is.null(blast_config)) {
    cli::cli_abort("No blast database configured in potato_config.yaml")
  }

  if (is.null(blast_config$files)) {
    cli::cli_abort("blast config must have 'files' array")
  }

  msg <- "Preparing BLAST annotation..."
  cli::cli_alert_info(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "info", message = msg)

  # Create filtered BLAST database from potato reference sequences
  blast_result <- create_blast_db(sack, potato_names)
  blast_db <- blast_result$blast_db
  all_subjects <- blast_result$all_subjects
  potatoes_with_blast <- blast_result$potatoes_with_blast
  fasta_content <- blast_result$fasta_content
  sack <- blast_result$sack

  # Get annotation session directory (created by create_blast_db)
  annotation_dir <- file.path(sack@sack_root, "results", "annotations", sack@metadata$annotation_session)

  # Build BLAST database
  msg <- "Building BLAST database..."
  cli::cli_alert_info(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "info", message = msg)
  makedb_cmd <- sprintf("makeblastdb -in %s -dbtype prot", shQuote(blast_db))
  if (!is.null(conda_env)) {
    makedb_cmd <- sprintf("conda run -n %s %s", conda_env, makedb_cmd)
  }
  system(makedb_cmd)
  msg <- "BLAST database ready"
  cli::cli_alert_success(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "success", message = msg)

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

  cli::cli_alert_info("Running BLAST on {length(genome_paths)} genome{?s}...")

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
    sprintf("%s run -n %s blastp -version 2>&1", conda_cmd, conda_env)
  } else {
    "blastp -version 2>&1"
  }
  tool_version <- tryCatch({
    version_output <- suppressWarnings(system(version_cmd, intern = TRUE))
    # Take first non-empty line
    version_line <- version_output[nchar(version_output) > 0][1]
    if (is.na(version_line)) "unknown" else version_line
  }, error = function(e) {
    "unknown"
  })

  # STEP 1: Run blast commands in parallel (just execute, return raw output + command)
  run_blast_cmd <- function(genome_path, genome_name, blast_db, conda_cmd, conda_env) {
    # BLAST tabular output format 6 (standard)
    # -evalue 1e-5: permissive threshold for performance (scoring phase can be stricter)
    # -max_target_seqs 500: limit hits per query for performance
    cmd <- sprintf("blastp -query %s -db %s -outfmt 6 -num_threads 1 -evalue 1e-5 -max_target_seqs 500",
                   shQuote(genome_path), shQuote(blast_db))

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
      format = "{cli::pb_spin} Running BLAST [{cli::pb_current}/{cli::pb_total}] | {cli::pb_bar} {cli::pb_percent}"
    ))
    progressr::with_progress({
      p <- progressr::progressor(along = genome_paths)
      results <- furrr::future_map(seq_along(genome_paths), function(i) {
        result <- run_blast_cmd(genome_paths[i], genome_names[i], blast_db, conda_cmd, conda_env)
        p()
        result
      }, .options = furrr::furrr_options(seed = TRUE))
    })
    future::plan(future::sequential)
  } else {
    cli::cli_progress_bar("Running BLAST", total = length(genome_paths))
    results <- purrr::map(seq_along(genome_paths), function(i) {
      cli::cli_progress_update()
      run_blast_cmd(genome_paths[i], genome_names[i], blast_db, conda_cmd, conda_env)
    })
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
    output_file <- file.path(annotation_dir, paste0(genome_names[i], ".blast.txt"))
    writeLines(raw_outputs[[i]], output_file)
  }

  # Write command log
  log_file <- file.path(annotation_dir, "blast.log")
  log_lines <- paste0(genome_names, "\t", commands)
  writeLines(c("genome\tcommand", log_lines), log_file)
  msg <- paste0("Saved outputs to ", annotation_dir)
  cli::cli_alert_success(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "success", message = msg)

  # STEP 2: Parse outputs sequentially using jakomics
  msg <- "Parsing BLAST results..."
  cli::cli_alert_info(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "info", message = msg)
  blast_module <- reticulate::import("jakomics.blast", delay_load = FALSE)

  blast_results <- purrr::map(raw_outputs, function(output_lines) {
    # Parse BLAST tabular output - returns ALL hits
    # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    parsed <- list()

    for (line in output_lines) {
      if (length(line) > 0) {
        fields <- strsplit(line, "\t")[[1]]
        if (length(fields) >= 12) {
          hit <- list(
            query = fields[1],
            subject = fields[2],
            pident = as.numeric(fields[3]),
            length = as.integer(fields[4]),
            mismatch = as.integer(fields[5]),
            gapopen = as.integer(fields[6]),
            qstart = as.integer(fields[7]),
            qend = as.integer(fields[8]),
            sstart = as.integer(fields[9]),
            send = as.integer(fields[10]),
            evalue = as.numeric(fields[11]),
            bitscore = as.numeric(fields[12])
          )

          # Group by subject ID
          if (is.null(parsed[[hit$subject]])) {
            parsed[[hit$subject]] <- list()
          }
          parsed[[hit$subject]][[length(parsed[[hit$subject]]) + 1]] <- hit
        }
      }
    }

    # Match to potato nodes
    blast_hits_to_tibble(parsed, potato_data, potato_hashes)
  })

  # Add to sack results
  sack@results$blast <- blast_results

  # Add provenance tracking
  # Build command templates with placeholders
  makedb_template <- sprintf("makeblastdb -in %s -dbtype prot", "{blast_db}")
  if (!is.null(conda_env)) {
    makedb_template <- sprintf("conda run -n %s %s", conda_env, makedb_template)
  }

  blastp_template <- "blastp -query {genome_path} -db {blast_db} -outfmt 6 -num_threads 1 -evalue 1e-5 -max_target_seqs 500"
  if (!is.null(conda_env)) {
    blastp_template <- sprintf("%s run -n %s %s", conda_cmd, conda_env, blastp_template)
  }

  # Store messages in provenance
  messages_tbl <- tibble::tibble(
    type = sapply(messages_list, function(x) x$type),
    message = sapply(messages_list, function(x) x$message)
  )

  sack@provenance$blast <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    tool_version = tool_version,
    conda_env = conda_env,
    workers = workers,
    potatoes_requested = names(potatoes),
    potatoes_with_genes = potatoes_with_blast,
    n_genomes = length(genome_paths),
    n_subjects = length(all_subjects),
    messages = messages_tbl,
    commands = list(
      fasta_content = fasta_content,
      blast_db = blast_db,
      makeblastdb_template = makedb_template,
      blastp_template = blastp_template,
      genome_paths = genome_paths
    )
  )

  msg <- "BLAST annotation complete"
  cli::cli_alert_success(msg)
  messages_list[[length(messages_list) + 1]] <- list(type = "success", message = msg)

  # Update stored messages
  sack@provenance$blast$messages <- tibble::tibble(
    type = sapply(messages_list, function(x) x$type),
    message = sapply(messages_list, function(x) x$message)
  )

  sack
}


#' Convert BLAST hits to tibble (internal)
#' @noRd
blast_hits_to_tibble <- function(parsed_hits, potato_data, potato_hashes) {
  # Build map of BLAST subject ID -> potato nodes
  subject_to_nodes <- list()
  for (potato_id in names(potato_data)) {
    potato <- potato_data[[potato_id]]
    for (node in potato$genes) {
      if (!is.null(node$databases$blast)) {
        for (subject_id in node$databases$blast) {
          if (is.null(subject_to_nodes[[subject_id]])) {
            subject_to_nodes[[subject_id]] <- list()
          }
          subject_to_nodes[[subject_id]][[length(subject_to_nodes[[subject_id]]) + 1]] <- list(
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
  for (subject_id in names(parsed_hits)) {
    hits_list <- parsed_hits[[subject_id]]
    nodes <- subject_to_nodes[[subject_id]]
    if (is.null(nodes)) next

    for (hit in hits_list) {
      for (node_info in nodes) {
        rows[[length(rows) + 1]] <- list(
          potato = node_info$potato,
          potato_hash = node_info$potato_hash,
          node_id = node_info$node_id,
          step = node_info$step,
          query = hit$query,
          subject = subject_id,
          pident = hit$pident,
          length = hit$length,
          evalue = hit$evalue,
          bitscore = hit$bitscore
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
      query = character(),
      subject = character(),
      pident = numeric(),
      length = integer(),
      evalue = numeric(),
      bitscore = numeric()
    ))
  }

  dplyr::bind_rows(rows)
}
