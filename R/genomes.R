#' Add genomes to a potato sack
#'
#' Register genome files with a potato sack without copying them. Can be called
#' multiple times to add genomes from different locations.
#'
#' @param sack PotatoSack object
#' @param path Character. Path to protein FASTA (.faa or .fasta) files. Can be:
#'   - A directory: adds all .faa/.fasta files in directory
#'   - A wildcard pattern: "~/data/mags/*.faa"
#'   - A single file: "~/data/genome1.faa"
#'   - A vector of files
#' @param validate Logical. If TRUE, validates .faa files (default: TRUE)
#' @param recursive Logical. If path is a directory, search recursively (default: FALSE)
#'
#' @returns Modified PotatoSack object with genomes added
#' @export

add_genomes <- function(sack, path, validate = TRUE, recursive = FALSE) {

  # Check jakomics is available
  if (!exists("jakomics")) {
    cli::cli_abort(c(
      "jakomics not loaded",
      "i" = "Check conda environment {.envvar potato} is activated"
    ))
  }

  # Handle different input types
  if (length(path) == 1 && dir.exists(path)) {
    # It's a directory
    genome_dir <- path
    genome_files <- NULL
  } else if (length(path) == 1 && grepl("\\*", path)) {
    # It's a wildcard pattern - expand it
    expanded <- Sys.glob(path)
    if (length(expanded) == 0) {
      cli::cli_abort("No files match pattern: {.file {path}}")
    }
    genome_dir <- NULL
    genome_files <- expanded
  } else {
    # It's a file or vector of files
    genome_dir <- NULL
    genome_files <- path
  }

  # Convert inputs for jakomics
  genome_dir_str <- if (is.null(genome_dir)) "" else genome_dir
  genome_files_list <- if (is.null(genome_files)) list() else as.list(genome_files)

  # Discover genome files via jakomics
  file_objects <- jakomics$utilities$get_files(
    file_names = genome_files_list,
    directory = genome_dir_str,
    file_type = list("faa", "fasta")
  )

  if (length(file_objects) == 0) {
    cli::cli_abort("No genome files found with {.file .faa} or {.file .fasta} extension")
  }

  cli::cli_alert_info("Found {length(file_objects)} genome file{?s}")

  # Validate FAA files if requested
  if (validate) {
    cli::cli_progress_bar("Validating genomes", total = length(file_objects))
    for (file_obj in file_objects) {
      cli::cli_progress_update()
      tryCatch({
        jakomics$utilities$validate_fasta(file_obj$file_path)
      }, error = function(e) {
        cli::cli_abort(c(
          "Invalid FASTA file: {.file {file_obj$name}}",
          "x" = e$message
        ))
      })
    }
    cli::cli_progress_done()
  }

  # Convert jakomics FILE objects to GenomeFile S7 objects
  new_genomes <- lapply(file_objects, jakomics_to_genome_file)

  # Check for duplicates (by file_path)
  if (length(sack@genomes) > 0) {
    existing_paths <- sapply(sack@genomes, function(g) g@file_path)
    new_paths <- sapply(new_genomes, function(g) g@file_path)

    duplicates <- intersect(existing_paths, new_paths)
    if (length(duplicates) > 0) {
      cli::cli_warn(c(
        "Skipping {length(duplicates)} genome{?s} already in sack",
        "i" = "{.file {basename(duplicates)}}"
      ))
      # Remove duplicates from new_genomes
      new_genomes <- new_genomes[!new_paths %in% duplicates]
    }
  }

  if (length(new_genomes) == 0) {
    cli::cli_alert_info("No new genomes added")
    return(sack)
  }

  # Add to @genomes slot
  sack@genomes <- c(sack@genomes, new_genomes)

  # Create tibble for @results (just genome names)
  new_results <- tibble::tibble(
    genome = sapply(new_genomes, function(g) g@short_name)
  )

  # Initialize or append to sack@results
  if (is.null(sack@results)) {
    sack@results <- new_results
  } else {
    sack@results <- dplyr::bind_rows(sack@results, new_results)
  }

  # Add provenance tracking (append, since add_genomes can be called multiple times)
  genome_provenance <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_added = length(new_genomes),
    genome_names = sapply(new_genomes, function(g) g@short_name),
    genome_paths = sapply(new_genomes, function(g) g@file_path)
  )

  if (is.null(sack@provenance$genomes)) {
    sack@provenance$genomes <- list(genome_provenance)
  } else {
    sack@provenance$genomes[[length(sack@provenance$genomes) + 1]] <- genome_provenance
  }

  cli::cli_alert_success("Added {length(new_genomes)} genome{?s}. Total: {length(sack@genomes)}")

  sack
}



