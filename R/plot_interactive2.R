#' Plot potato network with reaction-centric view (interactive visNetwork)
#'
#' Creates an interactive visNetwork visualization where nodes represent reactions.
#' Nodes sharing the same gene (KO ID) are connected with red dashed lines.
#'
#' @param potato Potato object or path to potato JSON file
#' @param sack PotatoSack object (optional, for genome detection status)
#' @param genome_name Genome name to show detection status (requires sack)
#' @param show_compounds Logical. Show compound nodes in bipartite graph (default: FALSE)
#' @param pathway For multi-pathway networks, show only this pathway (pathway ID). NULL = show all
#' @param layout Layout algorithm. "curated" uses x,y coordinates from potato JSON (default), or specify igraph layout: "fr", "kk", "dh", "circle", "nicely"
#' @param show_gene_connectors Logical. Show red dashed lines between nodes with same KO ID (default: TRUE)
#'
#' @return A visNetwork htmlwidget
#' @export
plot_potato_interactive2 <- function(potato, sack = NULL, genome_name = NULL,
                                      show_compounds = FALSE, pathway = NULL,
                                      layout = "curated", show_gene_connectors = TRUE) {

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg visNetwork} is required for interactive plots",
      "i" = "Install with: install.packages('visNetwork')"
    ))
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg igraph} is required for plotting")
  }

  # Prepare data
  prep <- prepare_potato_for_plotting(potato, sack, genome_name, show_compounds, pathway)

  # Calculate layout (scale up coordinates for visNetwork pixel space)
  node_coords <- calculate_node_layout(prep$potato, prep$g, prep$is_multi_pathway, show_compounds,
                                        layout = layout, scale_for_visnetwork = TRUE)

  # Build visNetwork with gene connectors
  build_visnetwork_with_gene_connectors(prep$potato, prep$g, node_coords, prep$node_status,
                                         prep$has_genome, show_gene_connectors)
}


#' Build visNetwork with gene connector lines (internal)
#' @noRd
build_visnetwork_with_gene_connectors <- function(potato, g, node_coords, node_status,
                                                   has_genome, show_gene_connectors = TRUE) {

  # Get edge list
  edge_list <- igraph::as_edgelist(g)
  edges_df <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2],
    arrows = "to",
    width = 2,
    dashes = FALSE,
    edge_type = "pathway",
    title = NA_character_,
    stringsAsFactors = FALSE
  )
  # Add color list column (replicate for each row)
  edges_df$color <- replicate(nrow(edges_df), list(color = "#2B7CE9", highlight = "#2B7CE9"), simplify = FALSE)

  # Add gene connector edges if requested
  if (show_gene_connectors) {
    # Get enzyme nodes only
    enzyme_nodes <- node_status[!node_status$is_compound_node, ]

    # Build map of node -> {ko_ids, step, gene_id}
    node_info <- list()
    for (i in seq_len(nrow(enzyme_nodes))) {
      node_name <- enzyme_nodes$name[i]
      gene_id <- enzyme_nodes$gene_id[i]

      # Get KO IDs
      gene_obj <- purrr::keep(potato@genes, ~ .x$id == gene_id)
      ko_ids <- character()
      if (length(gene_obj) > 0 && !is.null(gene_obj[[1]]$databases$kofam)) {
        ko_ids <- gene_obj[[1]]$databases$kofam
      }

      # Get step number from pathway
      step_num <- NA
      if (!is.null(potato@edges) && is.list(potato@edges)) {
        # Multi-pathway network
        for (pathway_id in names(potato@edges)) {
          pathway <- potato@edges[[pathway_id]]
          if (gene_id %in% names(pathway$nodes)) {
            step_num <- pathway$nodes[[gene_id]]$step
            break
          }
        }
      }

      node_info[[node_name]] <- list(
        gene_id = gene_id,
        ko_ids = ko_ids,
        step = step_num
      )
    }

    connector_edges <- list()

    # RED DASHED: Same gene (same KO), different reactions
    ko_map <- list()
    for (node_name in names(node_info)) {
      info <- node_info[[node_name]]
      for (ko in info$ko_ids) {
        if (is.null(ko_map[[ko]])) {
          ko_map[[ko]] <- character()
        }
        ko_map[[ko]] <- c(ko_map[[ko]], node_name)
      }
    }

    for (ko in names(ko_map)) {
      nodes_with_ko <- ko_map[[ko]]
      if (length(nodes_with_ko) > 1) {
        # Check if they're actually different gene IDs (avoid connecting OR branches)
        gene_ids <- sapply(nodes_with_ko, function(n) node_info[[n]]$gene_id)
        if (length(unique(gene_ids)) == 1) {
          # Same gene, multiple reactions - RED
          for (i in 1:(length(nodes_with_ko) - 1)) {
            for (j in (i + 1):length(nodes_with_ko)) {
              edge_df <- data.frame(
                from = nodes_with_ko[i],
                to = nodes_with_ko[j],
                arrows = "",
                width = 1,
                dashes = TRUE,
                edge_type = "same_gene",
                title = paste0("Same gene (", ko, "): different reactions"),
                stringsAsFactors = FALSE
              )
              edge_df$color <- list(list(color = "#FF0000", highlight = "#FF0000"))
              connector_edges[[length(connector_edges) + 1]] <- edge_df
            }
          }
        }
      }
    }

    # GREEN DASHED: Different genes, same step (OR branches)
    step_map <- list()
    for (node_name in names(node_info)) {
      info <- node_info[[node_name]]
      if (!is.na(info$step) && !is.list(info$step)) {
        step_key <- as.character(info$step)
        if (is.null(step_map[[step_key]])) {
          step_map[[step_key]] <- list()
        }
        step_map[[step_key]][[node_name]] <- info$gene_id
      }
    }

    for (step_num in names(step_map)) {
      nodes_at_step <- step_map[[step_num]]
      if (length(nodes_at_step) > 1) {
        # Check if they're different genes
        gene_ids <- unlist(nodes_at_step)
        if (length(unique(gene_ids)) > 1) {
          # Different genes, same step - GREEN
          node_names <- names(nodes_at_step)
          for (i in 1:(length(node_names) - 1)) {
            for (j in (i + 1):length(node_names)) {
              edge_df <- data.frame(
                from = node_names[i],
                to = node_names[j],
                arrows = "",
                width = 1,
                dashes = TRUE,
                edge_type = "or_branch",
                title = paste0("OR branch (step ", step_num, "): alternative enzymes"),
                stringsAsFactors = FALSE
              )
              edge_df$color <- list(list(color = "#00FF00", highlight = "#00FF00"))
              connector_edges[[length(connector_edges) + 1]] <- edge_df
            }
          }
        }
      }
    }

    if (length(connector_edges) > 0) {
      # Combine connector edges
      connector_dfs <- lapply(connector_edges, function(edge_df) {
        edge_df
      })
      connector_df <- do.call(rbind, connector_dfs)

      edges_df <- rbind(edges_df, connector_df)

      # Count by type
      red_count <- sum(connector_df$edge_type == "same_gene")
      green_count <- sum(connector_df$edge_type == "or_branch")
      if (red_count > 0) {
        cli::cli_alert_info("Added {red_count} red dashed edge(s) (same gene, different reactions)")
      }
      if (green_count > 0) {
        cli::cli_alert_info("Added {green_count} green dashed edge(s) (OR branches)")
      }
    }
  }

  # Check if using curated coordinates (already in screen space, typically 100-1000 range)
  # Layout algorithms produce coordinates in range ~-5 to 5, so we scale those up
  using_curated_coords <- any(!is.na(node_coords$x)) &&
                          max(abs(node_coords$x), na.rm = TRUE) > 50

  # Prepare nodes dataframe
  nodes_df <- data.frame(
    id = node_status$name,
    label = node_status$label,
    x = if (using_curated_coords) node_coords$x[match(node_status$name, node_coords$name)] else node_coords$x[match(node_status$name, node_coords$name)] * 100,
    y = if (using_curated_coords) node_coords$y[match(node_status$name, node_coords$name)] else node_coords$y[match(node_status$name, node_coords$name)] * 100,
    title = node_status$hover_text,
    shape = node_status$node_shape,
    stringsAsFactors = FALSE
  )

  # Color nodes by detection status (match original function styling)
  if (has_genome) {
    nodes_df$color.background <- sapply(seq_len(nrow(nodes_df)), function(i) {
      if (node_status$is_compound_node[i]) {
        if (node_status$status[i] == "Input") {
          return("#7FD399")  # Green for inputs
        } else if (node_status$status[i] == "Output") {
          return("#F77370")  # Red for outputs
        } else {
          return("#B3B3B3")  # Gray for intermediates
        }
      } else if (node_status$status[i] == "Detected") {
        return("#4CAF50")
      } else if (node_status$status[i] == "Partial") {
        return("#FFA726")
      } else if (node_status$status[i] == "Not detected") {
        return("#F44336")
      } else {
        return("#2196F3")  # Unknown
      }
    })
  } else {
    nodes_df$color.background <- sapply(seq_len(nrow(nodes_df)), function(i) {
      if (node_status$is_compound_node[i]) {
        if (node_status$status[i] == "Input") {
          return("#7FD399")  # Green for inputs
        } else if (node_status$status[i] == "Output") {
          return("#F77370")  # Red for outputs
        } else {
          return("#B3B3B3")  # Gray for intermediates
        }
      } else {
        return("#BBDEFB")  # Light blue for genes
      }
    })
  }

  nodes_df$color.border <- "#2B7CE9"
  nodes_df$color.highlight.background <- nodes_df$color.background
  nodes_df$color.highlight.border <- "#2B7CE9"

  # Build visNetwork (match original styling)
  vis <- visNetwork::visNetwork(nodes_df, edges_df, height = "100vh", width = "100%") %>%
    visNetwork::visNodes(
      font = list(size = 16),
      borderWidth = 2,
      shadow = TRUE
    ) %>%
    visNetwork::visEdges(
      smooth = FALSE,  # Straight edges like original
      arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5))
    ) %>%
    visNetwork::visPhysics(
      enabled = FALSE
    ) %>%
    visNetwork::visInteraction(
      dragNodes = TRUE,
      dragView = TRUE,
      zoomView = TRUE,
      navigationButtons = TRUE
    ) %>%
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = FALSE
    ) %>%
    visNetwork::visExport(
      type = "png",
      name = "potato_network",
      label = "Export as PNG"
    )

  vis
}
