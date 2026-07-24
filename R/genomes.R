#' Load genome files as jakomics FILE objects
#'
#' Discovers genome files (.faa, .gbk, .gb, .gbff) and returns them as a list of
#' jakomics FILE objects for downstream processing.
#'
#' @param genome_dir Directory containing genome files
#' @param genome_files Character vector of specific genome file paths (alternative to genome_dir)
#' @param validate Logical. If TRUE, validates .faa files using BioPython (default: TRUE)
#'
#' @returns List of jakomics FILE objects (reticulate Python objects)
#'
#' @details
#' This function discovers genome files and validates .faa files using
#' `jakomics.utilities.validate_fasta()`. Genbank files (.gb, .gbk, .gbff) are
#' not converted here - use `convert_genbank_to_faa()` for that.
#'
#' The returned FILE objects are reticulate Python objects with properties:
#' - `$file_path` - Full path to file
#' - `$name` - Filename with extension
#' - `$short_name` - Filename without extension
#' - `$suffix` - File extension (e.g., ".faa")
#' - `$id` - Unique ID (UUID hex)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load all genome files from directory
#' genomes <- prepare_genomes(genome_dir = "genomes/")
#'
#' # Load specific files
#' genomes <- prepare_genomes(genome_files = c("genome1.faa", "genome2.gbk"))
#'
#' # Access properties
#' genomes[[1]]$name
#' genomes[[1]]$file_path
#' }
prepare_genomes <- function(genome_dir = NULL, genome_files = NULL, validate = TRUE) {

  # Check jakomics is available
  if (!exists("jakomics")) {
    stop("jakomics not loaded. Check conda environment 'potato' is activated.", call. = FALSE)
  }

  # Input validation
  if (is.null(genome_dir) && is.null(genome_files)) {
    stop("Must provide either genome_dir or genome_files", call. = FALSE)
  }

  # Convert inputs
  genome_dir_str <- if (is.null(genome_dir)) "" else genome_dir
  genome_files_list <- if (is.null(genome_files)) list() else as.list(genome_files)

  # Discover genome files via jakomics
  file_objects <- jakomics$utilities$get_files(
    file_names = genome_files_list,
    directory = genome_dir_str,
    file_type = list("faa", "fasta", "gb", "gbk", "gbff")
  )

  if (length(file_objects) == 0) {
    stop("No genome files found matching extensions: faa, fasta, gb, gbk, gbff", call. = FALSE)
  }

  message("Found ", length(file_objects), " genome file(s)")

  # Validate FAA/FASTA files if requested
  if (validate) {
    for (file_obj in file_objects) {
      if (file_obj$suffix %in% c(".faa", ".fasta")) {
        message("  Validating: ", file_obj$name)
        tryCatch({
          jakomics$utilities$validate_fasta(file_obj$file_path)
        }, error = function(e) {
          stop("Invalid FASTA file: ", file_obj$name, "\n  ", e$message, call. = FALSE)
        })
      }
    }
  }

  file_objects
}


#' Convert genbank files to protein FASTA
#'
#' Takes jakomics FILE objects representing genbank files and converts them to
#' protein FASTA format using `jakomics.genome.GENOME.genbank_to_fasta()`.
#'
#' @param file_objects List of jakomics FILE objects (from `prepare_genomes()`)
#' @param output_dir Directory for output FAA files (default: "temp_genomes")
#' @param overwrite Logical. Overwrite existing FAA files (default: FALSE)
#'
#' @returns List of jakomics FILE objects pointing to the new FAA files
#'
#' @details
#' Only processes files with extensions .gb, .gbk, or .gbff. Other files are skipped
#' with a warning. Output files are named: `<genome_short_name>_<uuid>.faa`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load genomes (including genbank files)
#' genomes <- prepare_genomes(genome_dir = "genomes/")
#'
#' # Convert genbank to FAA
#' faa_genomes <- convert_genbank_to_faa(genomes)
#'
#' # Now all are FAA files
#' faa_genomes[[1]]$file_path
#' }
convert_genbank_to_faa <- function(file_objects, output_dir = "temp_genomes",
                                   overwrite = FALSE) {

  # Check jakomics is available
  if (!exists("jakomics")) {
    stop("jakomics not loaded. Check conda environment 'potato' is activated.", call. = FALSE)
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  faa_file_objects <- list()

  for (file_obj in file_objects) {
    # Skip non-genbank files
    if (!file_obj$suffix %in% c(".gb", ".gbk", ".gbff")) {
      warning("Skipping non-genbank file: ", file_obj$name)
      next
    }

    # Generate output path
    output_faa <- file.path(output_dir, paste0(file_obj$short_name, "_", file_obj$id, ".faa"))

    # Check if exists
    if (file.exists(output_faa) && !overwrite) {
      stop("Output file already exists: ", output_faa, "\n  Set overwrite=TRUE to replace", call. = FALSE)
    }

    message("Converting: ", file_obj$name, " -> ", basename(output_faa))

    # Convert via jakomics
    genome <- jakomics$genome$GENOME(file_obj)
    genome$genbank_to_fasta(write_faa = output_faa)

    # Create new FILE object for the FAA
    faa_file_obj <- jakomics$file$FILE(output_faa)
    faa_file_objects <- c(faa_file_objects, list(faa_file_obj))
  }

  if (length(faa_file_objects) == 0) {
    warning("No genbank files were converted")
  } else {
    message("Converted ", length(faa_file_objects), " genbank file(s)")
  }

  faa_file_objects
}
