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
    stop("jakomics not loaded. Check conda environment 'potato' is activated.", call. = FALSE)
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
      stop("No files match pattern: ", path, call. = FALSE)
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
    stop("No genome files found with .faa or .fasta extension (protein FASTA)", call. = FALSE)
  }

  message("Found ", length(file_objects), " genome file(s)")

  # Validate FAA files if requested
  if (validate) {
    for (file_obj in file_objects) {
      message("  Validating: ", file_obj$name)
      tryCatch({
        jakomics$utilities$validate_fasta(file_obj$file_path)
      }, error = function(e) {
        stop("Invalid FASTA file: ", file_obj$name, "\n  ", e$message, call. = FALSE)
      })
    }
  }

  new_genomes <- file_objects

  # Check for duplicates (by file_path)
  existing_paths <- sapply(sack@genomes, function(g) g$file_path)
  new_paths <- sapply(new_genomes, function(g) g$file_path)

  duplicates <- intersect(existing_paths, new_paths)
  if (length(duplicates) > 0) {
    warning("Skipping ", length(duplicates), " genome(s) already in sack:\n  ",
            paste(basename(duplicates), collapse = "\n  "), call. = FALSE)
    # Remove duplicates from new_genomes
    new_genomes <- new_genomes[!new_paths %in% duplicates]
  }

  if (length(new_genomes) == 0) {
    message("No new genomes added")
    return(sack)
  }

  # Append to existing genomes
  sack@genomes <- c(sack@genomes, new_genomes)

  message("Added ", length(new_genomes), " genome(s). Total: ", length(sack@genomes))

  sack
}



