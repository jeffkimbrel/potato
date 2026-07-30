#' Create step-based layout for potato graph (internal)
#' @noRd
create_step_layout <- function(potato, node_names) {

  # Extract step number from each node name (format: id_step)
  node_data <- purrr::map_dfr(node_names, function(node_name) {
    step <- as.integer(sub(".*_(\\d+)$", "\\1", node_name))
    node_id <- sub("_\\d+$", "", node_name)

    tibble::tibble(
      name = node_name,
      node_id = node_id,
      step = step
    )
  })

  # Group by step and assign x coordinates
  node_data <- node_data %>%
    dplyr::group_by(step) %>%
    dplyr::mutate(
      n_at_step = dplyr::n(),
      x_index = dplyr::row_number(),
      # Center alternatives horizontally: spread from -width/2 to +width/2
      x = ifelse(n_at_step == 1, 0, (x_index - 1) - (n_at_step - 1) / 2),
      # Y is negative step number (so pathway flows downward)
      y = -step
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(name, x, y)

  node_data
}


#' Plot potato pathway as a directed acyclic graph
#'
#' Visualizes pathway structure using ggraph, optionally highlighting detected
#' genes for a specific genome.
#'
#' @param potato Potato S7 object
#' @param sack Optional PotatoSack object (for highlighting detected genes)
#' @param genome_name Optional genome name (requires sack)
#' @param layout Layout algorithm: "sugiyama" (hierarchical), "fr" (force-directed), "tree"
#'
#' @export
plot_potato <- function(potato, sack = NULL, genome_name = NULL, layout = "sugiyama") {

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggraph} is required for plotting")
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg igraph} is required for plotting")
  }

  # Build igraph from potato edges
  g <- build_potato_graph(potato)

  # Get node detection status if genome provided
  if (!is.null(sack) && !is.null(genome_name)) {
    node_status <- get_node_status(potato, sack, genome_name)
  } else {
    node_status <- tibble::tibble(
      name = igraph::V(g)$name,
      detected = NA,
      status = "Unknown"
    )
  }

  # Add node labels (gene symbols without _step suffix)
  node_status <- node_status %>%
    dplyr::mutate(label = sub("_\\d+$", "", name))

  # Create manual layout based on step numbers
  # Y = -step (negative so it flows down), X = spread alternatives horizontally
  node_coords <- create_step_layout(potato, igraph::V(g)$name)

  # Create ggraph with manual layout
  p <- ggraph::ggraph(g, layout = "manual", x = node_coords$x, y = node_coords$y) +
    ggraph::geom_edge_link(
      arrow = grid::arrow(length = grid::unit(3, 'mm'), type = "closed"),
      end_cap = ggraph::circle(5, 'mm'),
      color = "gray50",
      alpha = 0.6
    ) +
    ggraph::geom_node_point(
      aes(color = node_status$status),
      size = 12
    ) +
    ggraph::geom_node_text(
      aes(label = node_status$label),
      size = 3,
      fontface = "bold"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Detected" = "#4CAF50",
        "Not detected" = "#F44336",
        "Unknown" = "#2196F3"
      ),
      na.value = "#2196F3"
    ) +
    ggplot2::labs(
      title = potato@name,
      subtitle = if (!is.null(genome_name)) paste("Genome:", genome_name) else NULL,
      color = "Detection Status"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20)) +
    potato_theme()

  p
}


#' Get edge labels (compounds) for plotting (internal)
#' @noRd
get_edge_labels <- function(potato, g) {

  # Get edge list from igraph
  edge_list <- igraph::as_edgelist(g)

  # Match edges to potato edges and extract compounds
  labels <- character(nrow(edge_list))

  for (i in seq_len(nrow(edge_list))) {
    from_node <- edge_list[i, 1]
    to_node <- edge_list[i, 2]

    # Find matching edge in potato
    matching_edge <- NULL
    for (edge in potato@edges) {
      if (edge$from == from_node && edge$to == to_node) {
        matching_edge <- edge
        break
      }
    }

    # Extract compound label
    if (!is.null(matching_edge) && !is.null(matching_edge$compound)) {
      labels[i] <- matching_edge$compound
    } else {
      labels[i] <- NA_character_
    }
  }

  labels
}


#' Get node detection status for plotting (internal)
#' @noRd
get_node_status <- function(potato, sack, genome_name) {

  # Find genome index
  genome_idx <- which(sack@results$genome == genome_name)
  if (length(genome_idx) == 0) {
    cli::cli_warn("Genome {.val {genome_name}} not found in results")
    return(tibble::tibble(
      name = igraph::V(build_potato_graph(potato))$name,
      detected = NA,
      status = "Unknown"
    ))
  }

  # Collect hits for this genome
  genome_hits <- list()

  if ("kofam" %in% names(sack@results)) {
    genome_hits$kofam <- sack@results$kofam[[genome_idx]]
  }
  if ("blast" %in% names(sack@results)) {
    genome_hits$blast <- sack@results$blast[[genome_idx]]
  }
  if ("hmm" %in% names(sack@results)) {
    genome_hits$hmm <- sack@results$hmm[[genome_idx]]
  }

  # Get all node names from graph
  g <- build_potato_graph(potato)
  node_names <- igraph::V(g)$name

  # Check each node
  status_data <- purrr::map_dfr(node_names, function(node_name) {
    # Extract node_id from node_name (remove _step suffix)
    node_id <- sub("_\\d+$", "", node_name)

    # Find corresponding node in potato
    node <- NULL
    for (n in potato@nodes) {
      if (n$id == node_id) {
        node <- n
        break
      }
    }

    if (is.null(node)) {
      return(tibble::tibble(
        name = node_name,
        detected = NA,
        status = "Unknown"
      ))
    }

    # Check if detected
    detected <- is_node_detected(node, potato@id, genome_hits)

    tibble::tibble(
      name = node_name,
      detected = detected,
      status = ifelse(detected, "Detected", "Not detected")
    )
  })

  status_data
}


#' Plot pathway completion heatmap across genomes
#'
#' Visualizes which pathways are present in which genomes as a heatmap using ggplot2.
#'
#' @param sack PotatoSack object with scores
#' @param cluster_rows Cluster pathways by similarity (default: TRUE)
#' @param cluster_cols Cluster genomes by similarity (default: TRUE)
#'
#' @export
plot_pathway_heatmap <- function(sack, cluster_rows = TRUE, cluster_cols = TRUE) {

  if (is.null(sack@scores)) {
    cli::cli_abort("No scores found. Run {.fn score_pathways} first.")
  }

  # Prepare data
  plot_data <- sack@scores %>%
    dplyr::select(genome, potato_name, present) %>%
    dplyr::mutate(present = ifelse(present, "Present", "Absent"))

  # Optionally cluster
  if (cluster_rows || cluster_cols) {
    # Create matrix for clustering
    mat <- sack@scores %>%
      dplyr::select(genome, potato_name, present) %>%
      tidyr::pivot_wider(names_from = potato_name, values_from = present) %>%
      tibble::column_to_rownames("genome") %>%
      as.matrix() * 1

    # Cluster rows (pathways)
    if (cluster_rows && nrow(mat) > 1) {
      row_order <- hclust(dist(t(mat)))$order
      pathway_levels <- colnames(mat)[row_order]
      plot_data$potato_name <- factor(plot_data$potato_name, levels = pathway_levels)
    }

    # Cluster cols (genomes)
    if (cluster_cols && ncol(mat) > 1) {
      col_order <- hclust(dist(mat))$order
      genome_levels <- rownames(mat)[col_order]
      plot_data$genome <- factor(plot_data$genome, levels = genome_levels)
    }
  }

  # Create heatmap
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = genome, y = potato_name, fill = present)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("Absent" = "gray90", "Present" = "#2E7D32"),
      name = "Status"
    ) +
    ggplot2::labs(
      title = "Pathway Presence/Absence Across Genomes",
      x = "Genome",
      y = "Pathway"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    potato_theme()

  p
}


#' Plot pathway completion fractions as barplot
#'
#' Shows completion fraction for each pathway in a specific genome using ggplot2.
#' Per-pathway thresholds are shown as markers on each bar.
#'
#' @param sack PotatoSack object with scores
#' @param genome_name Genome name to plot
#' @param show_thresholds Show per-pathway threshold markers (default: TRUE)
#'
#' @export
plot_genome_pathways <- function(sack, genome_name, show_thresholds = TRUE) {

  if (is.null(sack@scores)) {
    cli::cli_abort("No scores found. Run {.fn score_pathways} first.")
  }

  # Filter to this genome
  genome_scores <- sack@scores %>%
    dplyr::filter(genome == genome_name) %>%
    dplyr::arrange(dplyr::desc(fraction))

  if (nrow(genome_scores) == 0) {
    cli::cli_abort("No scores found for genome {.val {genome_name}}")
  }

  # Get threshold for each pathway
  genome_scores <- genome_scores %>%
    dplyr::mutate(
      threshold = purrr::map_dbl(potato, function(p_id) {
        potato <- sack@potatoes[[which(sapply(sack@potatoes, function(x) x@id == p_id))]]
        thresh <- potato@scoring$min_fraction
        if (is.null(thresh)) 0.75 else thresh
      })
    )

  # Order pathways by fraction
  genome_scores <- genome_scores %>%
    dplyr::mutate(potato_name = factor(potato_name, levels = potato_name))

  # Create plot
  p <- ggplot2::ggplot(genome_scores, ggplot2::aes(x = fraction, y = potato_name, fill = present)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#2E7D32", "FALSE" = "gray60"),
      labels = c("TRUE" = "Present", "FALSE" = "Absent"),
      name = "Status"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(
      title = paste("Pathway Completion in", genome_name),
      x = "Completion Fraction",
      y = NULL
    ) +
    potato_theme()

  # Add per-pathway threshold markers
  if (show_thresholds) {
    p <- p + ggplot2::geom_point(
      data = genome_scores,
      ggplot2::aes(x = threshold, y = potato_name),
      shape = 124,  # vertical bar
      size = 8,
      color = "red",
      alpha = 0.7,
      inherit.aes = FALSE
    )
  }

  p
}


#' Plot pathway completion summary
#'
#' Creates a summary plot showing number of pathways detected per genome.
#'
#' @param sack PotatoSack object with scores
#'
#' @export
plot_pathway_summary <- function(sack) {

  if (is.null(sack@scores)) {
    cli::cli_abort("No scores found. Run {.fn score_pathways} first.")
  }

  # Summarize pathways per genome
  summary_data <- sack@scores %>%
    dplyr::group_by(genome) %>%
    dplyr::summarize(
      total = dplyr::n(),
      present = sum(present),
      absent = total - present,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(present))

  # Reshape for stacking
  plot_data <- summary_data %>%
    tidyr::pivot_longer(cols = c(present, absent), names_to = "status", values_to = "count") %>%
    dplyr::mutate(
      genome = factor(genome, levels = summary_data$genome),
      status = factor(status, levels = c("present", "absent"))
    )

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = genome, y = count, fill = status)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("present" = "#2E7D32", "absent" = "gray60"),
      labels = c("present" = "Present", "absent" = "Absent"),
      name = "Status"
    ) +
    ggplot2::labs(
      title = "Pathways Detected Per Genome",
      x = "Genome",
      y = "Number of Pathways"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    potato_theme()

  p
}


#' Export potato graph to graphviz DOT format
#'
#' @param potato Potato S7 object
#' @param file Output file path
#'
#' @export
export_potato_dot <- function(potato, file) {

  g <- build_potato_graph(potato)

  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg igraph} is required")
  }

  # Write DOT file
  igraph::write_graph(g, file, format = "dot")

  cli::cli_alert_success("Wrote graphviz DOT to {.file {file}}")
}
