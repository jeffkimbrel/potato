#' Create filtered BLAST database from potato reference sequences
#'
#' Extracts only the reference sequences needed by potatoes from the configured
#' BLAST files, creates a filtered FASTA, and builds a BLAST database.
#'
#' @param sack PotatoSack object
#' @param potato_names Vector of potato names (NULL = all)
#'
#' @returns List with blast_db path and modified sack
#' @export
create_blast_db <- function(sack, potato_names = NULL) {

  # Filter potatoes
  if (is.null(potato_names)) {
    potatoes <- sack@potatoes
  } else {
    potatoes <- sack@potatoes[potato_names]
  }

  # Extract all unique BLAST subject IDs across potatoes
  all_subjects <- character()
  for (potato in potatoes) {
    # Handle both v1 (nodes) and v2 (genes)
    genes <- if (S7::S7_inherits(potato, PotatoV2)) {
      potato@genes
    } else {
      potato@genes
    }

    for (gene in genes) {
      if (!is.null(gene$databases$blast)) {
        all_subjects <- c(all_subjects, gene$databases$blast)
      }
    }
  }
  all_subjects <- unique(all_subjects)

  if (length(all_subjects) == 0) {
    cli::cli_abort("No blast terms found in selected potatoes")
  }

  cli::cli_alert_info("Creating BLAST database for {length(all_subjects)} reference sequence{?s}")

  # Get blast files from config
  blast_config <- sack@config$databases$blast
  blast_files <- blast_config$files

  if (is.null(blast_files) || length(blast_files) == 0) {
    cli::cli_abort("No blast files configured in potato_config.yaml")
  }

  # Read all FASTA files and extract needed sequences
  needed_seqs <- list()
  found_subjects <- character()

  for (fasta_file in blast_files) {
    if (!file.exists(fasta_file)) {
      cli::cli_warn("BLAST file not found: {.file {fasta_file}}")
      next
    }

    # Read FASTA file
    lines <- readLines(fasta_file)
    current_header <- NULL
    current_seq <- character()

    for (line in lines) {
      if (grepl("^>", line)) {
        # Save previous sequence if it was needed
        if (!is.null(current_header) && length(current_seq) > 0) {
          # Extract ID from header (first word after >)
          seq_id <- sub("^>([^ ]+).*", "\\1", current_header)
          if (seq_id %in% all_subjects) {
            needed_seqs[[seq_id]] <- list(
              header = current_header,
              sequence = paste(current_seq, collapse = "")
            )
            found_subjects <- c(found_subjects, seq_id)
          }
        }
        # Start new sequence
        current_header <- line
        current_seq <- character()
      } else {
        current_seq <- c(current_seq, line)
      }
    }

    # Don't forget the last sequence
    if (!is.null(current_header) && length(current_seq) > 0) {
      seq_id <- sub("^>([^ ]+).*", "\\1", current_header)
      if (seq_id %in% all_subjects) {
        needed_seqs[[seq_id]] <- list(
          header = current_header,
          sequence = paste(current_seq, collapse = "")
        )
        found_subjects <- c(found_subjects, seq_id)
      }
    }
  }

  # Check for missing sequences
  missing_subjects <- setdiff(all_subjects, found_subjects)
  if (length(missing_subjects) > 0) {
    cli::cli_warn(c(
      "Missing BLAST reference sequences:",
      "x" = "{.val {missing_subjects}}"
    ))
  }

  if (length(needed_seqs) == 0) {
    cli::cli_abort("No reference sequences found in BLAST files")
  }

  # Get or create annotation session timestamp
  if (is.null(sack@metadata$annotation_session)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    sack@metadata$annotation_session <- timestamp
  } else {
    timestamp <- sack@metadata$annotation_session
  }

  # Create results/annotations/{timestamp} directory if needed
  blast_dir <- file.path(sack@sack_root, "results", "annotations", timestamp)
  if (!dir.exists(blast_dir)) {
    dir.create(blast_dir, recursive = TRUE)
  }

  # Write filtered FASTA file
  filtered_fasta <- file.path(blast_dir, "filtered_blast_db.faa")
  fasta_lines <- character()
  for (seq_id in names(needed_seqs)) {
    fasta_lines <- c(fasta_lines,
                     needed_seqs[[seq_id]]$header,
                     needed_seqs[[seq_id]]$sequence)
  }
  writeLines(fasta_lines, filtered_fasta)

  cli::cli_alert_success("Created filtered FASTA with {length(needed_seqs)} sequence{?s}")

  list(blast_db = filtered_fasta, sack = sack)
}
