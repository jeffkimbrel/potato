#' Plot comprehensive pathway completion heatmap across all genomes
#'
#' Creates a heatmap showing pathway completion for all potatoes and genomes.
#' Tiles show fraction of pathway detected, optionally normalized by threshold.
#'
#' @param sack PotatoSack object with scores
#' @param normalize_by_threshold Logical. If TRUE, shows fraction/min_fraction ratio
#'   with color breakpoint at 1.0 (detection threshold). If FALSE, shows raw fraction.
#' @param cluster_rows Logical. Cluster pathways by similarity (default: FALSE)
#' @param cluster_cols Logical. Cluster genomes by similarity (default: TRUE)
#' @param show_labels Logical. Show pathway and genome labels (default: TRUE)
#' @param clustering_method Clustering method for hclust: "complete", "average", "single", "ward.D2" (default: "complete")
#'
#' @return A ggplot2 object
#' @export

plot_all_pathways_heatmap <- function(sack, normalize_by_threshold = TRUE,
                                       cluster_rows = FALSE, cluster_cols = TRUE,
                                       show_labels = TRUE, clustering_method = "complete") {


  if (is.null(sack@scores) || nrow(sack@scores) == 0) {
    cli::cli_abort("No scores found. Run {.fn score_pathways} first.")
  }

  # Prepare data
  scores_df <- sack@scores |>
    mutate(potato_count = n_distinct(potato), .by = "pathway_name") |> 
    dplyr::mutate(
      pathway_label = dplyr::case_when(
        potato_count == 1 ~ pathway_name,
    TRUE ~ paste0(potato, " : ", pathway_name)
      )
    ) |>
    select(-potato_count)

  # Calculate metric for heatmap
  if (normalize_by_threshold) {
    scores_df$value <- scores_df$fraction / scores_df$min_fraction
    value_label <- "Fraction / Threshold"
    breakpoint <- 1.0
  } else {
    scores_df$value <- scores_df$fraction
    value_label <- "Fraction Detected"
    breakpoint <- NULL
  }

  # Prepare matrix for optional clustering
  if (cluster_rows || cluster_cols) {
    # Create wide matrix (pathways × genomes)
    mat <- scores_df |>
      dplyr::select(pathway_label, genome, value) |>
      tidyr::pivot_wider(names_from = genome, values_from = value, values_fill = 0) |>
      tibble::column_to_rownames("pathway_label") |>
      as.matrix()

    # Cluster rows (pathways)
    if (cluster_rows) {
      row_clust <- hclust(dist(mat), method = clustering_method)
      row_order <- row_clust$order
      scores_df$pathway_label <- factor(scores_df$pathway_label,
                                         levels = rownames(mat)[row_order])
      cli::cli_alert_info("Clustered {nrow(mat)} pathways using {clustering_method} linkage")
    }

    # Cluster columns (genomes)
    if (cluster_cols) {
      col_clust <- hclust(dist(t(mat)), method = clustering_method)
      col_order <- col_clust$order
      scores_df$genome <- factor(scores_df$genome,
                                  levels = colnames(mat)[col_order])
      cli::cli_alert_info("Clustered {ncol(mat)} genomes using {clustering_method} linkage")
    }
  }

  # Create plot
  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = genome, y = pathway_label, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::labs(
      title = "Pathway Completion Across Genomes",
      x = "Genome",
      y = "Pathway",
      fill = value_label
    ) +
    potato_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      # axis.text.y = ggplot2::element_text(face = "plain"),
      # axis.title = ggplot2::element_text(face = "bold"),
      # plot.title = ggplot2::element_text(face = "bold"),
      # text = ggplot2::element_text(family = ""),  # Use default sans font
      panel.grid = ggplot2::element_blank()
    )

  # Color scale
  if (normalize_by_threshold) {
    # Diverging scale with midpoint at 1.0 (threshold)
    max_val <- max(scores_df$value, na.rm = TRUE)
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "#d73027",      # Red (below threshold)
        mid = "#fee090",      # Yellow (at threshold)
        high = "#1a9850",     # Green (above threshold)
        midpoint = 1.0,
        limits = c(0, max(max_val, 1.5)),
        breaks = c(0, 0.5, 1.0, 1.5),
        labels = c("0", "0.5", "1.0 (threshold)", "1.5"),
        na.value = "grey90"
      )
  } else {
    # Sequential scale for raw fraction
    p <- p +
      ggplot2::scale_fill_gradient(
        low = "#fee090",
        high = "#1a9850",
        limits = c(0, 1),
        na.value = "grey90"
      )
  }

  # Optional: hide labels for large heatmaps
  if (!show_labels) {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }

  p
}
