#' Load potato v2 schema
#' @param path Path to v2 potato JSON
#' @export
load_potato_v2 <- function(path) {
  potato_data <- jsonlite::read_json(path, simplifyVector = FALSE)

  # Validate schema version
  if (is.null(potato_data$schema_version) || potato_data$schema_version != "v2") {
    cli::cli_abort("Not a v2 schema potato (missing or wrong schema_version)")
  }

  cli::cli_alert_success("Loaded v2 potato: {potato_data$name}")
  cli::cli_alert_info("Genes: {length(potato_data$genes)}, Compounds: {length(potato_data$compounds)}")

  structure(potato_data, class = "potato_v2")
}


#' Build igraph from v2 potato
#' @param potato_v2 Potato loaded with load_potato_v2()
#' @param pathway_id Which pathway to build (default: first pathway)
#' @export
build_graph_v2 <- function(potato_v2, pathway_id = NULL) {

  if (is.null(pathway_id)) {
    pathway_id <- names(potato_v2$pathways)[1]
  }

  pathway <- potato_v2$pathways[[pathway_id]]
  cli::cli_alert_info("Building graph for pathway: {pathway$name}")

  # Get all edges
  edges <- pathway$edges

  # Collect edges by gene, then partition into distinct reactions
  gene_edges <- list()

  for (i in seq_along(edges)) {
    edge <- edges[[i]]
    from_node <- edge$from
    to_node <- edge$to

    # Edges where gene is the target (compound -> gene)
    if (grepl("^C\\d+", from_node) && !grepl("^C\\d+", to_node)) {
      gene_id <- to_node
      if (is.null(gene_edges[[gene_id]])) {
        gene_edges[[gene_id]] <- list(input_edges = list(), output_edges = list())
      }
      gene_edges[[gene_id]]$input_edges[[length(gene_edges[[gene_id]]$input_edges) + 1]] <- list(
        idx = i,
        compound = from_node
      )
    }

    # Edges where gene is the source (gene -> compound)
    if (!grepl("^C\\d+", from_node) && grepl("^C\\d+", to_node)) {
      gene_id <- from_node
      if (is.null(gene_edges[[gene_id]])) {
        gene_edges[[gene_id]] <- list(input_edges = list(), output_edges = list())
      }
      gene_edges[[gene_id]]$output_edges[[length(gene_edges[[gene_id]]$output_edges) + 1]] <- list(
        idx = i,
        compound = to_node
      )
    }
  }

  # Partition edges by (gene, reaction) pairs
  # Each unique gene+reaction combination becomes one node
  gene_instances <- list()
  edge_to_instance <- list()  # Map edge index to gene instance

  for (gene_id in names(gene_edges)) {
    input_edges <- gene_edges[[gene_id]]$input_edges
    output_edges <- gene_edges[[gene_id]]$output_edges

    # Group edges by reaction ID
    reaction_map <- list()

    for (e in input_edges) {
      rxn <- edges[[e$idx]]$reaction
      if (is.null(rxn)) rxn <- "default"
      if (is.null(reaction_map[[rxn]])) {
        reaction_map[[rxn]] <- list(inputs = list(), outputs = list())
      }
      reaction_map[[rxn]]$inputs[[length(reaction_map[[rxn]]$inputs) + 1]] <- e
    }

    for (e in output_edges) {
      rxn <- edges[[e$idx]]$reaction
      if (is.null(rxn)) rxn <- "default"
      if (is.null(reaction_map[[rxn]])) {
        reaction_map[[rxn]] <- list(inputs = list(), outputs = list())
      }
      reaction_map[[rxn]]$outputs[[length(reaction_map[[rxn]]$outputs) + 1]] <- e
    }

    # Create gene instance for each reaction
    if (length(reaction_map) == 1) {
      # Single reaction - use bare gene_id
      gene_instances[[gene_id]] <- list(
        gene_id = gene_id,
        instance_id = gene_id,
        reaction_id = names(reaction_map)[1],
        input_edges = input_edges,
        output_edges = output_edges
      )
      for (e in input_edges) edge_to_instance[[as.character(e$idx)]] <- gene_id
      for (e in output_edges) edge_to_instance[[as.character(e$idx)]] <- gene_id
    } else {
      # Multiple reactions - split by reaction
      for (rxn_id in names(reaction_map)) {
        rxn_data <- reaction_map[[rxn_id]]
        inst_id <- paste0(gene_id, "_", rxn_id)
        gene_instances[[inst_id]] <- list(
          gene_id = gene_id,
          instance_id = inst_id,
          reaction_id = rxn_id,
          input_edges = rxn_data$inputs,
          output_edges = rxn_data$outputs
        )
        for (e in rxn_data$inputs) edge_to_instance[[as.character(e$idx)]] <- inst_id
        for (e in rxn_data$outputs) edge_to_instance[[as.character(e$idx)]] <- inst_id
      }
      cli::cli_alert_info("Split {gene_id} into {length(reaction_map)} reaction nodes")
    }
  }

  # Build edge list using gene instances
  edge_list <- list()
  for (i in seq_along(edges)) {
    edge <- edges[[i]]
    from_node <- edge$from
    to_node <- edge$to

    # Map to gene instance if this edge involves a gene
    idx_str <- as.character(i)
    if (!grepl("^C\\d+", from_node) && !is.null(edge_to_instance[[idx_str]])) {
      from_node <- edge_to_instance[[idx_str]]
    }
    if (!grepl("^C\\d+", to_node) && !is.null(edge_to_instance[[idx_str]])) {
      to_node <- edge_to_instance[[idx_str]]
    }

    edge_list[[length(edge_list) + 1]] <- c(from_node, to_node)
  }

  # Create igraph
  edge_matrix <- do.call(rbind, edge_list)
  g <- igraph::graph_from_edgelist(edge_matrix, directed = TRUE)

  # Store original edges as graph attribute for reaction filtering
  g$pathway_edges <- edges

  # Add node attributes
  for (v in igraph::V(g)$name) {
    # Check if compound or gene
    if (grepl("^C\\d+", v)) {
      igraph::V(g)[v]$type <- "compound"
      # Find compound info
      compound <- purrr::keep(potato_v2$compounds, ~ .x$id == v)
      if (length(compound) > 0) {
        igraph::V(g)[v]$label <- compound[[1]]$name
      }
    } else {
      igraph::V(g)[v]$type <- "gene"
      # Strip _ctx or _R##### suffix to get original gene ID
      gene_id <- sub("_ctx\\d+$", "", v)
      gene_id <- sub("_R\\d+$", "", gene_id)
      gene <- purrr::keep(potato_v2$genes, ~ .x$id == gene_id)
      if (length(gene) > 0) {
        igraph::V(g)[v]$label <- gene[[1]]$name
        igraph::V(g)[v]$gene_id <- gene_id

        # Store databases, EC, reactions for hover text
        if (!is.null(gene[[1]]$databases)) {
          igraph::V(g)[v]$databases <- list(gene[[1]]$databases)
        }
        if (!is.null(gene[[1]]$ec)) {
          igraph::V(g)[v]$ec <- paste(gene[[1]]$ec, collapse = ", ")
        }
        if (!is.null(gene[[1]]$reactions)) {
          igraph::V(g)[v]$reactions <- paste(gene[[1]]$reactions, collapse = ", ")
        }
      }
    }
  }

  cli::cli_alert_success("Built graph: {igraph::vcount(g)} nodes, {igraph::ecount(g)} edges")

  g
}


#' Simple plot of v2 graph
#' @param potato_or_graph Either a potato_v2 object, path to v2 JSON, or igraph from build_graph_v2()
#' @param interactive Use interactive visNetwork plot (TRUE) or static ggraph (FALSE)
#' @export
plot_v2 <- function(potato_or_graph, interactive = TRUE) {

  # Build graph if needed
  if (is.character(potato_or_graph)) {
    # Path to JSON
    potato <- load_potato_v2(potato_or_graph)
    g <- build_graph_v2(potato)
  } else if (inherits(potato_or_graph, "potato_v2") ||
             (is.list(potato_or_graph) && !is.null(potato_or_graph$schema_version))) {
    # Potato object
    g <- build_graph_v2(potato_or_graph)
  } else {
    # Assume it's already a graph
    g <- potato_or_graph
  }

  if (interactive) {
    plot_v2_interactive(g)
  } else {
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggraph} required")
    }

    # Simple force-directed layout
    ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(3, "mm")),
                            end_cap = ggraph::circle(5, "mm"),
                            color = "gray50") +
      ggraph::geom_node_point(ggplot2::aes(color = type), size = 8) +
      ggraph::geom_node_text(ggplot2::aes(label = label), size = 3, repel = TRUE) +
      ggplot2::scale_color_manual(values = c("gene" = "#2196F3", "compound" = "#666666")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(title = "V2 Potato Graph")
  }
}


#' Interactive visNetwork plot of v2 graph
#' @param g Graph from build_graph_v2()
#' @export
plot_v2_interactive <- function(g) {
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg visNetwork} is required for interactive plots",
      "i" = "Install with: install.packages('visNetwork')"
    ))
  }

  # Get edge list
  edge_list <- igraph::as_edgelist(g)
  edges_df <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2],
    arrows = "to",
    width = 2,
    dashes = FALSE,
    stringsAsFactors = FALSE
  )
  edges_df$color <- replicate(nrow(edges_df), list(color = "#2B7CE9", highlight = "#2B7CE9"), simplify = FALSE)

  # Add red dashed edges between split gene nodes (same gene, different reactions)
  gene_nodes <- igraph::V(g)$name[igraph::V(g)$type == "gene"]
  gene_base_ids <- sub("_R\\d+$", "", gene_nodes)

  # Find genes with multiple reaction nodes
  for (base_id in unique(gene_base_ids)) {
    instances <- gene_nodes[gene_base_ids == base_id]
    if (length(instances) > 1) {
      # Connect all pairs of instances with red dashed line
      for (i in 1:(length(instances) - 1)) {
        for (j in (i + 1):length(instances)) {
          edges_df <- rbind(edges_df, data.frame(
            from = instances[i],
            to = instances[j],
            arrows = "none",
            width = 1,
            dashes = TRUE,
            color = I(list(list(color = "#FF0000", highlight = "#FF0000"))),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  # Add green dashed edges between alternative gene routes (OR branches)
  # Find genes that share same input compound(s) AND same output compound(s)
  for (gene1 in gene_nodes) {
    # Get predecessors and successors (compounds only)
    pred1 <- igraph::neighbors(g, gene1, mode = "in")
    pred1 <- pred1[igraph::V(g)[pred1]$type == "compound"]
    succ1 <- igraph::neighbors(g, gene1, mode = "out")
    succ1 <- succ1[igraph::V(g)[succ1]$type == "compound"]

    # Skip if no predecessors or successors
    if (length(pred1) == 0 || length(succ1) == 0) next

    # Check other genes
    for (gene2 in gene_nodes) {
      if (gene1 >= gene2) next  # Avoid duplicates and self

      pred2 <- igraph::neighbors(g, gene2, mode = "in")
      pred2 <- pred2[igraph::V(g)[pred2]$type == "compound"]
      succ2 <- igraph::neighbors(g, gene2, mode = "out")
      succ2 <- succ2[igraph::V(g)[succ2]$type == "compound"]

      # Check if they share same inputs and outputs
      if (length(pred2) > 0 && length(succ2) > 0 &&
          setequal(igraph::V(g)[pred1]$name, igraph::V(g)[pred2]$name) &&
          setequal(igraph::V(g)[succ1]$name, igraph::V(g)[succ2]$name)) {
        # These are alternative routes - add green dashed line
        edges_df <- rbind(edges_df, data.frame(
          from = gene1,
          to = gene2,
          arrows = "none",
          width = 1,
          dashes = TRUE,
          color = I(list(list(color = "#00AA00", highlight = "#00AA00"))),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Prepare nodes dataframe with hover text
  nodes_df <- data.frame(
    id = igraph::V(g)$name,
    label = igraph::V(g)$label,
    title = character(igraph::vcount(g)),
    shape = character(igraph::vcount(g)),
    stringsAsFactors = FALSE
  )

  # Build hover text
  for (i in seq_len(igraph::vcount(g))) {
    node_name <- nodes_df$id[i]
    node_type <- igraph::V(g)[node_name]$type

    if (node_type == "compound") {
      nodes_df$shape[i] <- "triangle"
      nodes_df$title[i] <- paste0(
        "Compound: ", igraph::V(g)[node_name]$label, "<br>",
        "ID: ", node_name
      )
    } else {
      nodes_df$shape[i] <- "dot"
      gene_id <- igraph::V(g)[node_name]$gene_id

      # Build hover text with all available info
      hover_parts <- c(
        paste0("Gene: ", gene_id),
        paste0("Enzyme: ", igraph::V(g)[node_name]$label)
      )

      # Add databases
      dbs <- igraph::V(g)[node_name]$databases
      if (!is.null(dbs) && length(dbs[[1]]) > 0) {
        db_lines <- sapply(names(dbs[[1]]), function(db_name) {
          terms <- paste(dbs[[1]][[db_name]], collapse = ", ")
          paste0(db_name, ": ", terms)
        })
        hover_parts <- c(hover_parts, db_lines)
      }

      # Add EC
      ec <- igraph::V(g)[node_name]$ec
      if (!is.null(ec) && !is.na(ec) && nchar(ec) > 0) {
        hover_parts <- c(hover_parts, paste0("EC: ", ec))
      }

      # Add reactions (filtered to pathway-relevant only)
      # Check if node name has _R##### suffix (split gene node)
      if (grepl("_R\\d+$", node_name)) {
        # Extract reaction from node name
        rxn_id <- sub(".*_(R\\d+)$", "\\1", node_name)
        hover_parts <- c(hover_parts, paste0("Reactions: ", rxn_id))
      } else {
        # Get reactions from edges involving this node
        pathway_edges <- g$pathway_edges
        relevant_rxns <- character()

        if (!is.null(pathway_edges)) {
          for (edge in pathway_edges) {
            # Check if this edge involves the current gene node
            if ((edge$from == node_name || edge$to == node_name) && !is.null(edge$reaction)) {
              relevant_rxns <- c(relevant_rxns, edge$reaction)
            }
          }
        }

        # If we found pathway-specific reactions, show only those
        if (length(relevant_rxns) > 0) {
          relevant_rxns <- unique(relevant_rxns)
          hover_parts <- c(hover_parts, paste0("Reactions: ", paste(relevant_rxns, collapse = ", ")))
        } else {
          # Fall back to all reactions from gene definition
          rxns <- igraph::V(g)[node_name]$reactions
          if (!is.null(rxns) && !is.na(rxns) && nchar(rxns) > 0) {
            hover_parts <- c(hover_parts, paste0("Reactions: ", rxns))
          }
        }
      }

      nodes_df$title[i] <- paste(hover_parts, collapse = "<br>")
    }
  }

  # Color nodes by type
  nodes_df$color.background <- ifelse(
    igraph::V(g)$type == "compound",
    "#B3B3B3",  # Gray for compounds
    "#BBDEFB"   # Light blue for genes
  )
  nodes_df$color.border <- "#2B7CE9"
  nodes_df$color.highlight.background <- nodes_df$color.background
  nodes_df$color.highlight.border <- "#2B7CE9"

  # Build visNetwork
  vis <- visNetwork::visNetwork(nodes_df, edges_df, height = "100vh", width = "100%") %>%
    visNetwork::visNodes(
      font = list(size = 16),
      borderWidth = 2,
      shadow = FALSE
    ) %>%
    visNetwork::visEdges(
      smooth = FALSE,
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
    )

  vis
}
