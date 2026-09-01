# Shared color palette for v2 plots
V2_COLORS <- list(
  gene = "#2196F3",        # Blue for genes
  compound = "#999999",    # Gray for compounds
  border = "#1976D2"       # Darker blue for borders
)

#' PotatoV2 S7 class
#'
#' Represents a v2 schema potato with genes, compounds, and pathways
#'
#' @export
PotatoV2 <- S7::new_class(
  "PotatoV2",
  properties = list(
    schema_version = S7::class_character,
    id = S7::class_character,
    name = S7::class_character,
    source = S7::class_character,
    tags = S7::class_character,
    notes = S7::class_character,
    genes = S7::class_list,
    compounds = S7::class_list,
    pathways = S7::class_list,
    json_path = S7::class_character
  ),
  validator = function(self) {
    if (self@schema_version != "v2") {
      "schema_version must be 'v2'"
    }
  }
)

#' Load potato v2 schema
#' @param path Path to v2 potato JSON
#' @export

load_potato <- function(path) {
  potato_data <- jsonlite::read_json(path, simplifyVector = FALSE)

  # Validate schema version
  if (is.null(potato_data$schema_version) || potato_data$schema_version != "v2") {
    cli::cli_abort("Not a v2 schema potato (missing or wrong schema_version)")
  }

  # Validate potato structure
  validation <- validate_potato(potato_data, strict = FALSE)
  if (!validation$valid) {
    cli::cli_abort(c(
      "Potato validation failed: {basename(path)}",
      "x" = "Errors found:",
      stats::setNames(validation$errors, rep("*", length(validation$errors)))
    ))
  }

  # Show warnings if any
  if (length(validation$warnings) > 0) {
    cli::cli_alert_warning("Potato has validation warnings:")
    for (warning in validation$warnings) {
      cli::cli_text("  {.emph {warning}}")
    }
  }

  cli::cli_alert_success("Loaded v2 potato: {potato_data$name}")
  cli::cli_alert_info("Genes: {length(potato_data$genes %||% list())}, Compounds: {length(potato_data$compounds %||% list())}")

  # Create S7 PotatoV2 object
  PotatoV2(
    schema_version = potato_data$schema_version,
    id = potato_data$id,
    name = potato_data$name,
    source = potato_data$source %||% "",
    tags = unlist(potato_data$tags) %||% character(0),
    notes = potato_data$notes %||% "",
    genes = potato_data$genes %||% list(),
    compounds = potato_data$compounds %||% list(),
    pathways = potato_data$pathways %||% list(),
    json_path = normalizePath(path, mustWork = FALSE)
  )
}


#' Build igraph from v2 potato
#' @param potato_v2 Potato loaded with load_potato()
#' @param pathway_id Which pathway to build (default: first pathway)
#' @export
build_graph <- function(potato_v2, pathway_id = NULL) {

  if (is.null(pathway_id)) {
    pathway_id <- names(potato_v2@pathways)[1]
  }

  pathway <- potato_v2@pathways[[pathway_id]]
  cli::cli_alert_info("Building graph for pathway: {pathway$name}")

  # Get all edges
  edges <- pathway$edges

  # Helper to check if node is a compound (not a gene)
  compound_ids <- sapply(potato_v2@compounds, function(c) c$id)
  is_compound <- function(node_id) {
    node_id %in% compound_ids
  }

  # Collect edges by gene, then partition into distinct reactions
  gene_edges <- list()

  for (i in seq_along(edges)) {
    edge <- edges[[i]]
    from_node <- edge$from
    to_node <- edge$to

    # Edges where gene is the target (compound -> gene)
    if (is_compound(from_node) && !is_compound(to_node)) {
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
    if (!is_compound(from_node) && is_compound(to_node)) {
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
    if (!is_compound(from_node) && !is.null(edge_to_instance[[idx_str]])) {
      from_node <- edge_to_instance[[idx_str]]
    }
    if (!is_compound(to_node) && !is.null(edge_to_instance[[idx_str]])) {
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
    if (is_compound(v)) {
      igraph::V(g)[v]$type <- "compound"
      # Find compound info
      compound <- purrr::keep(potato_v2@compounds, ~ .x$id == v)
      if (length(compound) > 0) {
        igraph::V(g)[v]$label <- compound[[1]]$name
      }
    } else {
      igraph::V(g)[v]$type <- "gene"
      # Strip _ctx or _R##### suffix to get original gene ID
      gene_id <- sub("_ctx\\d+$", "", v)
      gene_id <- sub("_R\\d+$", "", gene_id)
      gene <- purrr::keep(potato_v2@genes, ~ .x$id == gene_id)
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
#' @param potato_or_graph Either a potato_v2 object, path to v2 JSON, or igraph from build_graph()
#' @param interactive Use interactive visNetwork plot (TRUE) or static ggraph (FALSE)
#' @param layout Layout algorithm for static plot (e.g., "fr", "kk", "circle", "tree", "grid"). Default: "fr"
#' @export
plot_potato <- function(potato_or_graph, interactive = TRUE, layout = "fr") {

  # Build graph if needed
  potato <- NULL
  if (is.character(potato_or_graph)) {
    # Path to JSON
    potato <- load_potato(potato_or_graph)
    g <- build_graph(potato)
  } else if (inherits(potato_or_graph, "potato_v2") ||
             (is.list(potato_or_graph) && !is.null(potato_or_graph$schema_version))) {
    # Potato object
    potato <- potato_or_graph
    g <- build_graph(potato_or_graph)
  } else {
    # Assume it's already a graph
    g <- potato_or_graph
  }

  if (interactive) {
    plot_potato_interactive(g, layout = layout)
  } else {
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggraph} required")
    }

    # Get pathway name for title
    pathway_name <- if (!is.null(potato)) {
      # Get first pathway name
      first_pathway_id <- names(potato$pathways)[1]
      potato$pathways[[first_pathway_id]]$name
    } else {
      NULL
    }

    # Detect additional edges to add
    gene_nodes <- igraph::V(g)$name[igraph::V(g)$type == "gene"]
    gene_base_ids <- sub("_R\\d+$", "", gene_nodes)

    # Red dashed edges: split gene nodes (same gene, different reactions)
    red_edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
    for (base_id in unique(gene_base_ids)) {
      instances <- gene_nodes[gene_base_ids == base_id]
      if (length(instances) > 1) {
        for (i in 1:(length(instances) - 1)) {
          for (j in (i + 1):length(instances)) {
            red_edges <- rbind(red_edges, data.frame(from = instances[i], to = instances[j], stringsAsFactors = FALSE))
          }
        }
      }
    }

    # Green dashed edges: alternative routes (OR branches)
    green_edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
    for (gene1 in gene_nodes) {
      pred1 <- igraph::neighbors(g, gene1, mode = "in")
      pred1 <- pred1[igraph::V(g)[pred1]$type == "compound"]
      succ1 <- igraph::neighbors(g, gene1, mode = "out")
      succ1 <- succ1[igraph::V(g)[succ1]$type == "compound"]

      if (length(pred1) == 0 || length(succ1) == 0) next

      for (gene2 in gene_nodes) {
        if (gene1 >= gene2) next

        pred2 <- igraph::neighbors(g, gene2, mode = "in")
        pred2 <- pred2[igraph::V(g)[pred2]$type == "compound"]
        succ2 <- igraph::neighbors(g, gene2, mode = "out")
        succ2 <- succ2[igraph::V(g)[succ2]$type == "compound"]

        if (length(pred2) > 0 && length(succ2) > 0 &&
            setequal(igraph::V(g)[pred1]$name, igraph::V(g)[pred2]$name) &&
            setequal(igraph::V(g)[succ1]$name, igraph::V(g)[succ2]$name)) {
          green_edges <- rbind(green_edges, data.frame(from = gene1, to = gene2, stringsAsFactors = FALSE))
        }
      }
    }

    # Add additional edges to graph for plotting
    g_plot <- g
    # Always set edge_type attribute for all edges
    igraph::E(g_plot)$edge_type <- "pathway"

    if (nrow(red_edges) > 0) {
      g_plot <- g_plot + igraph::edges(as.vector(t(as.matrix(red_edges))), edge_type = "split_gene")
    }
    if (nrow(green_edges) > 0) {
      g_plot <- g_plot + igraph::edges(as.vector(t(as.matrix(green_edges))), edge_type = "alternative")
    }

    # Apply layout - pathway edges with arrows, dashed edges without
    p <- ggraph::ggraph(g_plot, layout = layout) +
      # Pathway edges (solid, with arrows)
      ggraph::geom_edge_link(ggplot2::aes(filter = edge_type == "pathway"),
                            color = "gray50",
                            arrow = grid::arrow(length = grid::unit(3, "mm")),
                            end_cap = ggraph::circle(5, "mm")) +
      # Red dashed edges (split genes, no arrows)
      ggraph::geom_edge_link(ggplot2::aes(filter = edge_type == "split_gene"),
                            color = "#FF0000",
                            linetype = "dashed",
                            arrow = NULL) +
      # Green dashed edges (alternatives, no arrows)
      ggraph::geom_edge_link(ggplot2::aes(filter = edge_type == "alternative"),
                            color = "#00AA00",
                            linetype = "dashed",
                            arrow = NULL) +
      ggraph::geom_node_point(ggplot2::aes(color = type, shape = type), size = 8) +
      ggraph::geom_node_text(ggplot2::aes(label = label), size = 3, repel = TRUE) +
      ggplot2::scale_color_manual(values = c("gene" = V2_COLORS$gene, "compound" = V2_COLORS$compound)) +
      ggplot2::scale_shape_manual(values = c("gene" = 19, "compound" = 17)) +  # 19 = circle, 17 = triangle
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.position = "none"
      )

    # Add title if we have pathway name
    if (!is.null(pathway_name)) {
      p <- p + ggplot2::labs(title = pathway_name)
    }

    p
  }
}


#' Interactive visNetwork plot of v2 graph
#' @param g Graph from build_graph()
#' @param layout Layout algorithm (e.g., "fr", "kk", "circle", "tree", "grid"). Default: "fr"
#' @param height Height of the plot (default: "100vh")
#' @export
plot_potato_interactive <- function(g, layout = "fr", height = "100vh") {
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
  edges_df$color <- replicate(nrow(edges_df), list(color = "#216bde", highlight = "#216bde"), simplify = FALSE)

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

  # Calculate layout positions
  layout_func <- switch(layout,
    "fr" = igraph::layout_with_fr,
    "kk" = igraph::layout_with_kk,
    "circle" = igraph::layout_in_circle,
    "tree" = igraph::layout_as_tree,
    "grid" = igraph::layout_on_grid,
    "sugiyama" = function(g) igraph::layout_with_sugiyama(g)$layout,
    igraph::layout_with_fr  # default to fr
  )
  coords <- layout_func(g)

  # Prepare nodes dataframe with hover text
  nodes_df <- data.frame(
    id = igraph::V(g)$name,
    label = igraph::V(g)$label,
    title = character(igraph::vcount(g)),
    shape = character(igraph::vcount(g)),
    x = coords[, 1] * 100,  # Scale for visNetwork
    y = coords[, 2] * 100,
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
              # Handle both single strings and arrays
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

  # Color nodes by type (using shared palette)
  nodes_df$color.background <- ifelse(
    igraph::V(g)$type == "compound",
    V2_COLORS$compound,  # Gray for compounds
    V2_COLORS$gene       # Blue for genes
  )
  nodes_df$color.border <- V2_COLORS$border
  nodes_df$color.highlight.background <- nodes_df$color.background
  nodes_df$color.highlight.border <- V2_COLORS$border

  # Build visNetwork
  vis <- visNetwork::visNetwork(nodes_df, edges_df, height = height, width = "100%") %>%
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
