#' Get detailed blast annotation results
#'
#' Returns detailed information about blast annotation hits, including
#' per-genome and per-potato summaries and plots.
#'
#' @param sack PotatoSack object
#'
#' @return A list containing:
#'   \item{summary}{Tibble with per-genome hit counts}
#'   \item{per_potato}{Tibble with per-potato hit statistics}
#'   \item{plot}{ggplot object showing hits per genome}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' details <- get_blast_details(sack)
#' details$summary
#' details$per_potato
#' details$plot
#' }
get_blast_details <- function(sack) {

  # Check if blast has been run
  if (is.null(sack@results$blast)) {
    cli::cli_alert_info("BLAST annotation has not been run on this sack")
    return(list(
      summary = tibble::tibble(),
      per_potato = tibble::tibble(),
      plot = NULL
    ))
  }

  # Extract blast results (list of tibbles, one per genome)
  blast_list <- sack@results$blast

  # Get genome names
  genome_names <- sapply(sack@genomes, function(g) g@short_name)

  # Build per-genome summary
  summary_rows <- list()
  for (i in seq_along(blast_list)) {
    hits <- blast_list[[i]]
    summary_rows[[i]] <- list(
      genome = genome_names[i],
      n_hits = nrow(hits),
      n_unique_subjects = length(unique(hits$subject)),
      n_potatoes = length(unique(hits$potato))
    )
  }

  summary_tbl <- tibble::tibble(
    genome = sapply(summary_rows, function(x) x$genome),
    n_hits = sapply(summary_rows, function(x) x$n_hits),
    n_unique_subjects = sapply(summary_rows, function(x) x$n_unique_subjects),
    n_potatoes = sapply(summary_rows, function(x) x$n_potatoes)
  )

  # Build per-potato summary (across all genomes)
  # Add genome names to each tibble
  blast_with_genomes <- list()
  for (i in seq_along(blast_list)) {
    hits <- blast_list[[i]]
    hits$genome <- genome_names[i]
    blast_with_genomes[[i]] <- hits
  }
  all_hits <- dplyr::bind_rows(blast_with_genomes)

  per_potato_tbl <- all_hits |>
    dplyr::group_by(potato) |>
    dplyr::summarise(
      n_hits = dplyr::n(),
      n_unique_subjects = dplyr::n_distinct(subject),
      n_genomes_with_hits = dplyr::n_distinct(genome),
      .groups = "drop"
    )

  # Create plot - hits per genome
  p <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = genome, y = n_hits)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Genome", y = "BLAST Hits") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(
    summary = summary_tbl,
    per_potato = per_potato_tbl,
    plot = p
  )
}
