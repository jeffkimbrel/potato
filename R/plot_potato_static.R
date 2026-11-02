#' Plot potato network with pathway hulls (static ggraph)
#'
#' Creates a static ggplot2/ggraph visualization of a multi-pathway potato network
#' with convex hulls around each pathway's genes.
#'
#' @param potato Potato V2 object or path to potato JSON file
#' @param pathway For multi-pathway networks, show only specific pathway(s) (default: NULL = all)
#' @param layout Layout algorithm: "fr", "kk", "sugiyama", "tree", "circle" (default: "fr")
#' @param show_compounds Show compound nodes (default: TRUE)
#'
#' @return A ggplot2 object
#' @export

plot_potato_static <- function(potato, pathway = NULL, layout = "fr", show_compounds = TRUE) {

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggraph} is required for static plots. Install with: install.packages('ggraph')")
  }

  if (!requireNamespace("ggforce", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggforce} is required for pathway hulls. Install with: install.packages('ggforce')")
  }

  # Load if path provided
  if (is.character(potato)) {
    potato <- load_potato(potato)
  }

  is_multi_pathway <- !is.null(potato@pathways) && length(potato@pathways) > 1

  # Determine which pathways to show
  if (is.null(pathway)) {
    # Show all pathways
    pathway_ids <- names(potato@pathways)
  } else {
    # Show specific pathway(s)
    pathway_ids <- pathway
    # Validate pathway IDs
    invalid <- pathway_ids[!pathway_ids %in% names(potato@pathways)]
    if (length(invalid) > 0) {
      cli::cli_abort("Pathway(s) not found: {paste(invalid, collapse = ', ')}")
    }
  }

  # Build combined graph from all requested pathways
  all_nodes <- character()
  all_edges <- list()

  for (pid in pathway_ids) {
    pw <- potato@pathways[[pid]]

    # Collect nodes from edges
    if (!is.null(pw$edges)) {
      for (edge in pw$edges) {
        all_nodes <- c(all_nodes, edge$from, edge$to)
        all_edges[[length(all_edges) + 1]] <- list(from = edge$from, to = edge$to)
      }
    }
  }

  all_nodes <- unique(all_nodes)

  # Filter to genes/compounds based on show_compounds
  compound_ids <- sapply(potato@compounds, function(c) c$id)
  gene_ids <- sapply(potato@genes, function(g) g$id)

  if (!show_compounds) {
    all_nodes <- all_nodes[all_nodes %in% gene_ids]
    # Filter edges to gene→gene only
    all_edges <- Filter(function(e) e$from %in% gene_ids && e$to %in% gene_ids, all_edges)
  }

  # Build igraph
  if (length(all_edges) == 0) {
    cli::cli_abort("No edges to plot")
  }

  edge_df <- data.frame(
    from = sapply(all_edges, function(e) e$from),
    to = sapply(all_edges, function(e) e$to)
  )

  # Create vertices dataframe with type
  vertices_df <- data.frame(
    name = all_nodes,
    type = ifelse(all_nodes %in% compound_ids, "compound", "gene"),
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(edge_df, directed = TRUE, vertices = vertices_df)

  # Calculate layout
  layout_func <- switch(layout,
    "fr" = igraph::layout_with_fr,
    "kk" = igraph::layout_with_kk,
    "circle" = igraph::layout_in_circle,
    "tree" = igraph::layout_as_tree,
    "sugiyama" = function(g) igraph::layout_with_sugiyama(g)$layout,
    igraph::layout_with_fr  # default
  )

  coords <- layout_func(g)
  node_layout <- data.frame(
    name = igraph::V(g)$name,
    x = coords[, 1],
    y = coords[, 2],
    type = igraph::V(g)$type
  )

  # Create base ggraph plot
  p <- ggraph::ggraph(g, layout = "manual", x = node_layout$x, y = node_layout$y)

  # Add pathway convex hulls if multi-pathway
  if (is_multi_pathway && length(pathway_ids) > 1) {
    pathway_hull_data <- list()

    # Generate colors for pathways
    pathway_colors <- grDevices::rainbow(length(pathway_ids), alpha = 0.3)
    names(pathway_colors) <- pathway_ids

    for (i in seq_along(pathway_ids)) {
      pid <- pathway_ids[i]
      pw <- potato@pathways[[pid]]
      pw_name <- pw$name %||% pid

      # Get nodes for this pathway
      pw_node_ids <- character()
      if (!is.null(pw$genes)) {
        # Explicit genes list
        pw_node_ids <- pw$genes
      }
      # Also add from edges
      if (!is.null(pw$edges)) {
        for (edge in pw$edges) {
          pw_node_ids <- c(pw_node_ids, edge$from, edge$to)
        }
      }
      pw_node_ids <- unique(pw_node_ids)

      # Filter to nodes in graph
      pw_node_ids <- pw_node_ids[pw_node_ids %in% node_layout$name]

      if (length(pw_node_ids) > 0) {
        pw_coords <- node_layout[node_layout$name %in% pw_node_ids, ]
        pw_coords$pathway_id <- pid
        pw_coords$pathway_name <- pw_name
        pw_coords$pathway_color <- pathway_colors[i]
        pathway_hull_data[[i]] <- pw_coords
      }
    }

    if (length(pathway_hull_data) > 0) {
      hull_df <- do.call(rbind, pathway_hull_data)

      p <- p +
        ggforce::geom_mark_hull(
          data = hull_df,
          ggplot2::aes(x = x, y = y, fill = pathway_name, color = pathway_name),
          alpha = 0.15,
          expand = grid::unit(8, "mm"),
          radius = grid::unit(10, "mm"),
          concavity = 0.5,
          show.legend = TRUE,
          size = 1.2
        ) +
        ggplot2::scale_fill_discrete(name = "Pathway") +
        ggplot2::scale_color_discrete(name = "Pathway")
    }
  }

  # Add edges
  p <- p +
    ggraph::geom_edge_link(
      arrow = grid::arrow(length = grid::unit(3, 'mm'), type = "closed"),
      end_cap = ggraph::circle(3, 'mm'),
      color = "gray50",
      alpha = 0.6,
      linewidth = 0.8
    )

  # Add nodes
  p <- p +
    ggraph::geom_node_point(
      ggplot2::aes(shape = type, color = type),
      size = 8,
      stroke = 1.5
    ) +
    ggplot2::scale_shape_manual(
      values = c("gene" = 19, "compound" = 17),  # circle, triangle
      name = "Type"
    ) +
    ggplot2::scale_color_manual(
      values = c("gene" = V2_COLORS$gene, "compound" = V2_COLORS$compound),
      name = "Type"
    )

  # Add labels
  p <- p +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      size = 3,
      repel = TRUE,
      max.overlaps = 20
    )

  # Theme
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::labs(title = potato@name)

  return(p)
}
