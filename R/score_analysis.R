#' Identify genes missing across genomes
#'
#' Analyzes which genes are systematically missing, suggesting issues with
#' database representation or thresholds.
#'
#' @param sack PotatoSack object with scores
#' @param potato_name Optional: filter to specific pathway
#' @param min_genomes Minimum genomes to consider (default: 3)
#'
#' @returns Tibble showing genes and how often they're missing
#' @export
summarize_missing_genes <- function(sack, potato_name = NULL, min_genomes = 3) {

  if (is.null(sack@results)) {
    cli::cli_abort("No annotation results found. Run annotation tools first.")
  }

  # Filter to specific potato if requested
  potatoes <- sack@potatoes
  if (!is.null(potato_name)) {
    potatoes <- potatoes[sapply(potatoes, function(p) p@name == potato_name)]
  }

  # Track missing genes
  missing_data <- list()

  for (genome_name in sack@results$genome) {
    genome_idx <- which(sack@results$genome == genome_name)

    # Collect hits
    genome_hits <- list()
    if ("kofam" %in% names(sack@results)) genome_hits$kofam <- sack@results$kofam[[genome_idx]]
    if ("blast" %in% names(sack@results)) genome_hits$blast <- sack@results$blast[[genome_idx]]
    if ("hmm" %in% names(sack@results)) genome_hits$hmm <- sack@results$hmm[[genome_idx]]

    # Check each potato
    for (potato in potatoes) {
      for (node in potato@nodes) {
        detected <- is_node_detected(node, potato@id, genome_hits)

        missing_data[[length(missing_data) + 1]] <- list(
          potato = potato@name,
          gene_id = node$id,
          gene_name = node$name,
          genome = genome_name,
          detected = detected
        )
      }
    }
  }

  # Convert to tibble and summarize
  missing_df <- dplyr::bind_rows(missing_data) %>%
    dplyr::group_by(potato, gene_id, gene_name) %>%
    dplyr::summarize(
      total_genomes = dplyr::n(),
      times_detected = sum(detected),
      times_missing = sum(!detected),
      fraction_missing = times_missing / total_genomes,
      .groups = "drop"
    ) %>%
    dplyr::filter(total_genomes >= min_genomes) %>%
    dplyr::arrange(dplyr::desc(fraction_missing))

  missing_df
}


#' Find pathways that are close to being called present
#'
#' Identifies "near miss" pathways that are just below the threshold,
#' which may indicate threshold tuning issues.
#'
#' @param sack PotatoSack object with scores
#' @param buffer Distance from threshold to consider "near" (default: 0.1)
#'
#' @returns Tibble of near-miss pathways
#' @export
find_near_miss_pathways <- function(sack, buffer = 0.1) {

  if (is.null(sack@scores)) {
    cli::cli_abort("No scores found. Run {.fn score_pathways} first.")
  }

  # Find pathways that are close but not present
  near_miss <- sack@scores %>%
    dplyr::filter(!present) %>%
    dplyr::mutate(
      # Get threshold for each pathway
      threshold = purrr::map_dbl(potato, function(p_id) {
        potato <- sack@potatoes[[which(sapply(sack@potatoes, function(x) x@id == p_id))]]
        thresh <- potato@scoring$min_fraction
        if (is.null(thresh)) 0.75 else thresh
      }),
      distance_from_threshold = threshold - fraction,
      near_miss = distance_from_threshold <= buffer
    ) %>%
    dplyr::filter(near_miss) %>%
    dplyr::select(genome, potato_name, fraction, threshold,
                  distance_from_threshold, steps_detected, steps_total) %>%
    dplyr::arrange(distance_from_threshold)

  if (nrow(near_miss) == 0) {
    cli::cli_alert_info("No near-miss pathways found (within {buffer} of threshold)")
  } else {
    cli::cli_alert_success("Found {nrow(near_miss)} near-miss pathway{?s}")
  }

  near_miss
}


#' Plot near-miss pathways
#'
#' Visualizes pathways that are close to being called present, highlighting
#' which are just below the threshold.
#'
#' @param sack PotatoSack object with scores
#' @param genome_name Optional: filter to specific genome
#' @param buffer Distance from threshold to highlight (default: 0.1)
#'
#' @export
plot_near_miss_pathways <- function(sack, genome_name = NULL, buffer = 0.1) {

  if (is.null(sack@scores)) {
    cli::cli_abort("No scores found. Run {.fn score_pathways} first.")
  }

  # Prepare data
  plot_data <- sack@scores %>%
    dplyr::mutate(
      # Get threshold for each pathway
      threshold = purrr::map_dbl(potato, function(p_id) {
        potato <- sack@potatoes[[which(sapply(sack@potatoes, function(x) x@id == p_id))]]
        thresh <- potato@scoring$min_fraction
        if (is.null(thresh)) 0.75 else thresh
      }),
      distance_from_threshold = threshold - fraction,
      status = dplyr::case_when(
        present ~ "Present",
        distance_from_threshold <= buffer ~ "Near miss",
        TRUE ~ "Absent"
      )
    )

  # Filter to genome if specified
  if (!is.null(genome_name)) {
    plot_data <- plot_data %>% dplyr::filter(genome == genome_name)
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fraction, y = potato_name,
                                                 color = status, shape = status)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(aes(xintercept = threshold),
                        linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(
      values = c("Present" = "#2E7D32", "Near miss" = "#FF9800", "Absent" = "gray60"),
      name = "Status"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Present" = 19, "Near miss" = 17, "Absent" = 1),
      name = "Status"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title = if (!is.null(genome_name)) paste("Pathway Completion:", genome_name) else "Pathway Completion",
      subtitle = paste0("Near-miss buffer: ", buffer),
      x = "Completion Fraction",
      y = "Pathway"
    ) +
    potato_theme()

  # Facet by genome if showing multiple
  if (is.null(genome_name) && length(unique(plot_data$genome)) > 1) {
    p <- p + ggplot2::facet_wrap(~genome, ncol = 3)
  }

  p
}
