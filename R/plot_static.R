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

  # For multi-pathway networks: assign pathway colors to nodes
  if (prep$is_multi_pathway) {
    # Generate colors for pathways
    pathway_names <- sapply(prep$potato@edges, function(p) p$name %||% "")
    pathway_colors <- jakR2::palette_jak(n = length(pathway_names), p = "sunset")
    names(pathway_colors) <- pathway_names

    # Assign each enzyme node to its pathway(s)
    prep$node_status$pathway_for_color <- sapply(seq_len(nrow(prep$node_status)), function(i) {
      node_name <- prep$node_status$name[i]

      # Compounds get gray (neutral)
      if (prep$node_status$is_compound_node[i]) {
        return("Compound")
      }

      # Find which pathway this gene belongs to
      for (pathway_id in names(prep$potato@edges)) {
        pathway <- prep$potato@edges[[pathway_id]]
        if (node_name %in% names(pathway$nodes)) {
          return(pathway$name %||% pathway_id)
        }
      }
      return("Unknown")
    })
  } else {
    # Single pathway: use detection status for enzymes, gray for compounds
    prep$node_status$pathway_for_color <- ifelse(
      prep$node_status$is_compound_node,
      "Compound",
      prep$node_status$status
    )
  }

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

      if (requireNamespace("ggforce", quietly = TRUE)) {
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
      # ggplot2::aes(color = prep$node_status$pathway_for_color, shape = prep$node_status$shape_code),
      ggplot2::aes(shape = prep$node_status$shape_code), color = "gray70",
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

  # Color scales - different logic for multi-pathway vs single-pathway
  if (prep$is_multi_pathway) {
    # Multi-pathway: nodes colored by pathway, compounds gray
    pathway_names <- sapply(prep$potato@edges, function(p) p$name %||% "")
    pathway_colors <- jakR2::palette_jak(n = length(pathway_names), p = "sunset")
    names(pathway_colors) <- pathway_names

    # Add gray for compounds
    all_colors <- c(pathway_colors, "Compound" = "#B3B3B3")

    p <- p +
      ggplot2::scale_color_manual(
        values = all_colors,
        name = "Pathway"
      ) +
      ggplot2::scale_fill_manual(
        values = pathway_colors,
        name = "Pathway"
      ) +
      ggplot2::labs(
        title = prep$potato@name,
        subtitle = if (prep$has_genome) paste("Genome:", genome_name) else NULL
      )
  } else {
    # Single-pathway: use detection status colors
    if (prep$has_genome) {
      p <- p +
        ggplot2::scale_color_manual(
          values = c(
            "Detected" = "#7FD399",
            "Partial" = "#FFB84D",
            "Not detected" = "#F77370",
            "Unknown" = "#6AB8F5",
            "Compound" = "#B3B3B3"
          ),
          na.value = "#6AB8F5",
          name = "Detection Status"
        ) +
        ggplot2::labs(
          title = prep$potato@name,
          subtitle = paste("Genome:", genome_name)
        )
    } else {
      p <- p +
        ggplot2::scale_color_manual(
          values = c(
            "Unknown" = "#6AB8F5",
            "Compound" = "#B3B3B3"
          ),
          na.value = "#6AB8F5",
          guide = "none"
        ) +
        ggplot2::labs(
          title = prep$potato@name,
          subtitle = NULL
        )
    }
  }

  p
}
