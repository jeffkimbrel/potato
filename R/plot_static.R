#' Plot potato network (static ggraph)
#'
#' Creates a static ggplot2/ggraph visualization of a potato pathway network.
#' Supports curated layouts, pathway filtering, and genome detection status.
#'
#' @param potato Potato object or path to potato JSON file
#' @param sack PotatoSack object (optional, for genome detection status)
#' @param genome_name Genome name to show detection status (requires sack)
#' @param show_compounds Logical. Show compound nodes in bipartite graph (default: FALSE)
#' @param layout Layout algorithm: "fr", "kk", "sugiyama", "tree", "circle", "grid" (default: "fr")
#' @param pathway For multi-pathway networks, show only this pathway (pathway ID). NULL = show all
#'
#' @return A ggplot2 object
#' @export

plot_potato_static <- function(potato, sack = NULL, genome_name = NULL,
                                show_compounds = FALSE, layout = "fr", pathway = NULL) {

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggraph} is required for static plots")
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg igraph} is required for plotting")
  }

  # Prepare data
  prep <- prepare_potato_for_plotting(potato, sack, genome_name, show_compounds, pathway)

  # Calculate layout
  node_coords <- calculate_node_layout(prep$potato, prep$g, prep$is_multi_pathway, show_compounds, layout)

  # Map shapes to numeric values for ggraph (solid shapes - no border)
  # 19 = solid circle, 17 = solid triangle, 18 = solid diamond, 15 = solid square
  shape_map <- c("circle" = 19, "triangle" = 17, "diamond" = 18, "square" = 15)
  prep$node_status$shape_code <- shape_map[prep$node_status$node_shape]

  # Assign color status to each node
  # Nodes are ONLY colored by detection status (never by pathway)
  # Pathway colors are reserved for hulls only
  prep$node_status$color_status <- ifelse(
    prep$node_status$is_compound_node,
    "Compound",  # Compounds always gray
    if (prep$has_genome) prep$node_status$status else "NoGenome"  # Detection status or generic gray
  )

  # Create ggraph with manual layout - flip Y axis to match visNetwork screen coordinates
  p <- ggraph::ggraph(prep$g, layout = "manual", x = node_coords$x, y = -node_coords$y)

  # Add pathway convex hulls for multi-pathway networks
  if (prep$is_multi_pathway) {
    pathway_nodes <- list()

    # Use the already-generated pathway colors
    pathway_names <- sapply(prep$potato@edges, function(p) p$name %||% "")
    pathway_colors <- jakR2::palette_jak(n = length(pathway_names), p = "sunset")
    names(pathway_colors) <- pathway_names

    for (pathway_id in names(prep$potato@edges)) {
      pathway <- prep$potato@edges[[pathway_id]]
      node_ids <- names(pathway$nodes)
      pathway_name <- pathway$name %||% pathway_id
      pathway_color <- pathway_colors[pathway_name]

      for (node_id in node_ids) {
        idx <- which(node_coords$name == node_id)
        if (length(idx) > 0) {
          pathway_nodes[[length(pathway_nodes) + 1]] <- data.frame(
            x = node_coords$x[idx],
            y = -node_coords$y[idx],  # Flip Y to match ggraph coordinates
            pathway = pathway_id,
            pathway_name = pathway_name,
            pathway_color = pathway_color
          )
        }
      }
    }

    if (length(pathway_nodes) > 0) {
      pathway_hull_data <- do.call(rbind, pathway_nodes)

      # Filter to pathways with 2+ genes (geom_mark_hull needs multiple points)
      pathway_gene_counts <- table(pathway_hull_data$pathway_name)
      pathway_hull_data <- pathway_hull_data[pathway_hull_data$pathway_name %in% names(pathway_gene_counts[pathway_gene_counts >= 2]), ]

      if (nrow(pathway_hull_data) > 0 && requireNamespace("ggforce", quietly = TRUE)) {
        p <- p +
          ggforce::geom_mark_hull(
            data = pathway_hull_data,
            ggplot2::aes(x = x, y = y, fill = pathway_name, color = pathway_name),
            alpha = 0.1,
            expand = grid::unit(12, "mm"),
            radius = grid::unit(13, "mm"),
            concavity = 0.5,
            show.legend = TRUE,
            size = 1.5
          )
      } else if (nrow(pathway_hull_data) == 0) {
        cli::cli_alert_info("All pathways are single-gene; no hulls to draw (nodes colored by pathway)")
      } else {
        cli::cli_warn("Package {.pkg ggforce} is required for pathway hulls")
      }
    }
  }

  # Add edges
  p <- p +
    ggraph::geom_edge_link(
      arrow = grid::arrow(length = grid::unit(2, 'mm'), type = "closed"),
      end_cap = ggraph::circle(2.25, 'mm'),  # Match node size of 6 (was 3mm for size 8)
      color = "gray50",
      alpha = 0.6
    )

  # Add nodes (solid shapes - no separate border)
  p <- p +
    ggraph::geom_node_point(
      ggplot2::aes(color = prep$node_status$color_status, shape = prep$node_status$shape_code),
      size = 6
    )

  # Add text labels with ggrepel
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = data.frame(
          x = node_coords$x,
          y = -node_coords$y,
          label = prep$node_status$label
        ),
        ggplot2::aes(x = x, y = y, label = label),
        size = 4,
        family = "sans",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "gray70",
        segment.size = 0.3,
        max.overlaps = Inf
      )
  } else {
    p <- p +
      ggraph::geom_node_text(
        ggplot2::aes(label = prep$node_status$label),
        size = 4,
        family = "sans"
      )
  }

  # Add scales and theme
  p <- p +
    ggplot2::scale_shape_identity() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20)) +
    potato_theme()

  # Build color palette for nodes (detection status only, never pathway)
  node_colors <- c(
    "Detected" = "#7FD399",
    "Partial" = "#FFB84D",
    "Not detected" = "#F77370",
    "Unknown" = "#6AB8F5",
    "NoGenome" = "#B3B3B3",   # Gray for enzymes when no genome
    "Compound" = "#B3B3B3"     # Gray for compounds
  )

  # Apply node color scale
  if (prep$has_genome) {
    p <- p +
      ggplot2::scale_color_manual(
        values = node_colors,
        breaks = c("Detected", "Partial", "Not detected", "Unknown", "Compound"),
        name = "Detection Status",
        na.value = "#6AB8F5"
      )
  } else {
    # No genome: all nodes gray, no legend needed
    p <- p +
      ggplot2::scale_color_manual(
        values = node_colors,
        guide = "none"
      )
  }

  # Add title/subtitle
  p <- p +
    ggplot2::labs(
      title = prep$potato@name,
      subtitle = if (prep$has_genome) paste("Genome:", genome_name) else NULL
    )

  # Add pathway hull fill scale for multi-pathway networks
  if (prep$is_multi_pathway) {
    pathway_names <- sapply(prep$potato@edges, function(p) p$name %||% "")
    pathway_colors <- jakR2::palette_jak(n = length(pathway_names), p = "sunset")
    names(pathway_colors) <- pathway_names

    p <- p +
      ggplot2::scale_fill_manual(
        values = pathway_colors,
        name = "Pathway"
      )
  }

  p
}
