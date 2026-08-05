#' Plot pathway prevalence across genomes
#'
#' Creates a horizontal bar chart showing how many genomes have each pathway.
#' Useful for identifying widespread vs rare metabolic capabilities.
#'
#' @param sack PotatoSack object with scoring results
#' @param min_genomes Minimum number of genomes to include pathway (default: 0)
#' @param pathway_type For multi-pathway networks, filter to "variant" or "independent" (default: NULL = all)
#' @param sort_by Sort bars by "count" (default), "name", or "fraction"
#'
#' @return A ggplot2 object
#' @export
plot_pathway_prevalence <- function(sack, min_genomes = 0, pathway_type = NULL, sort_by = "count") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required")
  }

  if (is.null(sack@scores) || nrow(sack@scores) == 0) {
    cli::cli_abort("No scoring results found. Run score_pathways() first.")
  }

  # Calculate pathway prevalence
  total_genomes <- length(unique(sack@scores$genome))

  # For multi-pathway networks, use pathway_name, otherwise use potato_name
  if ("pathway_name" %in% names(sack@scores)) {
    # Multi-pathway: count by pathway_name
    prevalence <- sack@scores %>%
      dplyr::filter(present == TRUE) %>%
      dplyr::group_by(potato, potato_name, pathway_name) %>%
      dplyr::summarise(
        genomes_detected = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        fraction_genomes = genomes_detected / total_genomes,
        label = pathway_name
      )

    # Add pathway type if available
    if (!is.null(sack@potatoes)) {
      pathway_types <- list()
      for (pot_name in names(sack@potatoes)) {
        pot <- sack@potatoes[[pot_name]]
        if (!is.null(pot@edges) && is.list(pot@edges) && length(pot@edges) > 0) {
          for (pathway_id in names(pot@edges)) {
            pathway <- pot@edges[[pathway_id]]
            pathway_types[[pathway$name %||% pathway_id]] <- pathway$type %||% "unknown"
          }
        }
      }
      prevalence$type <- sapply(prevalence$pathway_name, function(n) pathway_types[[n]] %||% "single")
    }
  } else {
    # Single-pathway: count by potato_name
    prevalence <- sack@scores %>%
      dplyr::filter(present == TRUE) %>%
      dplyr::group_by(potato, potato_name) %>%
      dplyr::summarise(
        genomes_detected = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        fraction_genomes = genomes_detected / total_genomes,
        label = potato_name,
        type = "single"
      )
  }

  # Filter by pathway type if specified
  if (!is.null(pathway_type) && "type" %in% names(prevalence)) {
    prevalence <- prevalence %>%
      dplyr::filter(type == pathway_type)
  }

  # Filter by minimum genomes
  prevalence <- prevalence %>%
    dplyr::filter(genomes_detected >= min_genomes)

  if (nrow(prevalence) == 0) {
    cli::cli_warn("No pathways meet the filtering criteria")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # Sort
  if (sort_by == "count") {
    prevalence <- prevalence %>%
      dplyr::arrange(genomes_detected)
  } else if (sort_by == "fraction") {
    prevalence <- prevalence %>%
      dplyr::arrange(fraction_genomes)
  } else if (sort_by == "name") {
    prevalence <- prevalence %>%
      dplyr::arrange(label)
  }

  prevalence$label <- factor(prevalence$label, levels = prevalence$label)

  # Create plot
  p <- ggplot2::ggplot(prevalence, ggplot2::aes(x = genomes_detected, y = label)) +
    ggplot2::geom_col(ggplot2::aes(fill = type), width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%d (%.0f%%)", genomes_detected, fraction_genomes * 100)),
      hjust = -0.1,
      size = 3
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15)),
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::labs(
      title = "Pathway Prevalence Across Genomes",
      subtitle = sprintf("%d genomes analyzed", total_genomes),
      x = "Number of genomes with pathway",
      y = NULL,
      fill = "Type"
    ) +
    potato_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = if ("type" %in% names(prevalence) && length(unique(prevalence$type)) > 1) "right" else "none"
    )

  # Color scale
  if ("type" %in% names(prevalence)) {
    type_colors <- c(
      "single" = "#7570b3",
      "variant" = "#e7298a",
      "independent" = "#66a61e",
      "unknown" = "#999999"
    )
    p <- p + ggplot2::scale_fill_manual(values = type_colors)
  }

  p
}


#' Plot pathway uniqueness distribution
#'
#' Creates a histogram showing the distribution of pathway prevalence.
#' Shows how many pathways are found in 1 genome, 2 genomes, etc.
#'
#' @param sack PotatoSack object with scoring results
#'
#' @return A ggplot2 object
#' @export
plot_pathway_uniqueness <- function(sack) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required")
  }

  if (is.null(sack@scores) || nrow(sack@scores) == 0) {
    cli::cli_abort("No scoring results found. Run score_pathways() first.")
  }

  # Calculate pathway prevalence
  if ("pathway_name" %in% names(sack@scores)) {
    # Multi-pathway networks
    prevalence <- sack@scores %>%
      dplyr::filter(present == TRUE) %>%
      dplyr::group_by(pathway_name) %>%
      dplyr::summarise(genomes_detected = dplyr::n(), .groups = "drop")
  } else {
    # Single-pathway potatoes
    prevalence <- sack@scores %>%
      dplyr::filter(present == TRUE) %>%
      dplyr::group_by(potato_name) %>%
      dplyr::summarise(genomes_detected = dplyr::n(), .groups = "drop")
  }

  # Count pathways by genome count
  distribution <- prevalence %>%
    dplyr::group_by(genomes_detected) %>%
    dplyr::summarise(pathway_count = dplyr::n(), .groups = "drop")

  # Categorize
  total_genomes <- length(unique(sack@scores$genome))
  distribution <- distribution %>%
    dplyr::mutate(
      category = dplyr::case_when(
        genomes_detected == 1 ~ "Unique (1 genome)",
        genomes_detected <= total_genomes * 0.25 ~ "Rare (< 25%)",
        genomes_detected <= total_genomes * 0.75 ~ "Common (25-75%)",
        TRUE ~ "Widespread (> 75%)"
      )
    )

  # Create plot
  p <- ggplot2::ggplot(distribution, ggplot2::aes(x = genomes_detected, y = pathway_count, fill = category)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = pathway_count),
      vjust = -0.5,
      size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Unique (1 genome)" = "#e41a1c",
        "Rare (< 25%)" = "#ff7f00",
        "Common (25-75%)" = "#4daf4a",
        "Widespread (> 75%)" = "#377eb8"
      )
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      title = "Pathway Uniqueness Distribution",
      subtitle = sprintf("%d pathways across %d genomes", nrow(prevalence), total_genomes),
      x = "Number of genomes",
      y = "Number of pathways",
      fill = "Category"
    ) +
    potato_theme() +
    ggplot2::theme(legend.position = "right")

  p
}
