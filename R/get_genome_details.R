#' Get detailed genome information
#'
#' Returns a list with detailed information about all genomes in the sack,
#' including file sizes, protein counts, and QC plots.
#'
#' @param sack PotatoSack object
#'
#' @return A list containing:
#'   \item{summary}{Tibble with genome information (genome, file_path, file_type, md5, file_size_mb, n_proteins, added_in_call)}
#'   \item{plot}{ggplot object with protein count bar plot}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- get_genome_details(sack)
#' result$summary
#' result$plot
#' }

get_genome_details <- function(sack) {

  if (length(sack@genomes) == 0) {
    cli::cli_alert_info("No genomes in sack")
    empty_summary <- tibble::tibble(
      genome = character(),
      file_path = character(),
      file_type = character(),
      md5 = character(),
      file_size_mb = numeric(),
      n_proteins = integer(),
      added_in_call = integer()
    )
    return(list(
      summary = empty_summary,
      plot = NULL
    ))
  }

  # Build details tibble
  details <- tibble::tibble(
    genome = sapply(sack@genomes, function(g) g@short_name),
    file_path = sapply(sack@genomes, function(g) g@file_path),
    file_type = sapply(sack@genomes, function(g) g@file_type),
    md5 = sapply(sack@genomes, function(g) g@md5)
  )

  # Add file sizes (MB)
  details$file_size_mb <- sapply(details$file_path, function(path) {
    if (file.exists(path)) {
      file.info(path)$size / 1024^2  # Convert bytes to MB
    } else {
      NA_real_
    }
  })

  # Count proteins (number of ">" lines in FASTA)
  details$n_proteins <- sapply(details$file_path, function(path) {
    if (file.exists(path)) {
      # Count lines starting with ">"
      length(grep("^>", readLines(path, warn = FALSE)))
    } else {
      NA_integer_
    }
  })

  # Add provenance information if available
  if (!is.null(sack@provenance$genomes)) {
    # Build mapping of genome name -> call number
    genome_call_map <- list()

    for (call_num in seq_along(sack@provenance$genomes)) {
      prov <- sack@provenance$genomes[[call_num]]
      for (genome_name in prov$genome_names) {
        # First occurrence wins (in case genome was re-added, which should be prevented)
        if (is.null(genome_call_map[[genome_name]])) {
          genome_call_map[[genome_name]] <- call_num
        }
      }
    }

    # Add to details
    details$added_in_call <- sapply(details$genome, function(g) {
      call_num <- genome_call_map[[g]]
      if (is.null(call_num)) NA_integer_ else call_num
    })
  } else {
    details$added_in_call <- NA_integer_
  }

  # Create plot - protein count only
  p <- ggplot2::ggplot(details, ggplot2::aes(x = genome, y = n_proteins)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Genome", y = "Protein Count") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(
    summary = details,
    plot = p
  )
}
