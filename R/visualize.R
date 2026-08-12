#' Normalize compound names by sorting multiple compounds
#' @noRd
normalize_compound_name <- function(compound_name) {
  if (is.null(compound_name)) {
    return(compound_name)
  }
  if (is.na(compound_name) || nchar(compound_name) == 0) {
    return(compound_name)
  }

  # Split by " + " (with spaces), sort, rejoin
  parts <- strsplit(compound_name, " \\+ ")[[1]]
  if (length(parts) > 1) {
    parts <- sort(trimws(parts))
    return(paste(parts, collapse = " + "))
  }
  return(compound_name)
}

#' Build bipartite graph with compound nodes (internal)
#' @noRd
build_bipartite_graph <- function(potato) {

  # Check if this is a multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  # Collect all unique compounds from edges
  compound_info <- list()
  compound_status <- list()  # Track if compound is input, output, or intermediate

  # Collect all edges (with proper DAG node IDs for networks)
  all_dag_edges <- list()

  if (is_network) {
    # Multi-pathway network - use gene IDs directly (not step-based)
    for (pathway_id in names(potato@edges)) {
      pathway <- potato@edges[[pathway_id]]

      # Process edges for this pathway - use gene IDs directly
      if (!is.null(pathway$edges) && length(pathway$edges) > 0) {
        for (edge in pathway$edges) {
        from_id <- edge$from
        to_id <- edge$to

        if (!is.null(edge$compound)) {
          # Parse compound string (may contain multiple compounds)
          # Normalize compound name by sorting parts
          compound_name <- normalize_compound_name(edge$compound)
          kegg_id <- edge$kegg_compound

          # Use KEGG ID for deduplication if available
          if (!is.null(kegg_id) && !is.na(kegg_id) && nchar(kegg_id) > 0) {
            compound_id <- paste0("COMPOUND_", kegg_id)
          } else {
            compound_id <- paste0("COMPOUND_", gsub("[^A-Za-z0-9]", "_", compound_name))
          }

          if (!(compound_id %in% names(compound_info))) {
            compound_info[[compound_id]] <- list(
              id = compound_id,
              name = compound_name,
              kegg_id = kegg_id
            )
            compound_status[[compound_id]] <- "intermediate"
          }

          # Only add edges if both endpoints exist (skip external compounds with null endpoints)
          if (!is.null(from_id)) {
            all_dag_edges[[length(all_dag_edges) + 1]] <- list(
              from = from_id,
              to = compound_id,
              pathway = pathway_id,
              pathway_name = pathway$name %||% pathway_id,
              pathway_type = pathway$type
            )
          }
          if (!is.null(to_id)) {
            all_dag_edges[[length(all_dag_edges) + 1]] <- list(
              from = compound_id,
              to = to_id,
              pathway = pathway_id,
              pathway_name = pathway$name %||% pathway_id,
              pathway_type = pathway$type
            )
          }
        } else {
          # Direct edge - gene to gene (skip if either endpoint is null)
          if (!is.null(from_id) && !is.null(to_id)) {
            all_dag_edges[[length(all_dag_edges) + 1]] <- list(
              from = from_id,
              to = to_id,
              pathway = pathway_id,
              pathway_name = pathway$name %||% pathway_id,
              pathway_type = pathway$type
            )
          }
        }
        }
      }

      # Add input compound for this pathway (don't split, use full string)
      if (!is.null(pathway$input)) {
        compound_name <- pathway$input$compound
        kegg_id <- pathway$input$kegg_compound

        # Use KEGG ID for deduplication if available
        if (!is.null(kegg_id) && !is.na(kegg_id) && nchar(kegg_id) > 0) {
          compound_id <- paste0("COMPOUND_", kegg_id)
        } else {
          compound_id <- paste0("COMPOUND_", gsub("[^A-Za-z0-9]", "_", compound_name))
        }

        # Check if already exists as intermediate - keep as input (green priority)
        if (!(compound_id %in% names(compound_info))) {
          compound_info[[compound_id]] <- list(
            id = compound_id,
            name = compound_name,
            kegg_id = kegg_id
          )
          compound_status[[compound_id]] <- "input"
        } else if (compound_status[[compound_id]] == "intermediate") {
          # Upgrade from intermediate to input
          compound_status[[compound_id]] <- "input"
        }

        # Connect input compound to target genes
        if (!is.null(pathway$input$targets)) {
          for (target in pathway$input$targets) {
            all_dag_edges[[length(all_dag_edges) + 1]] <- list(
              from = compound_id,
              to = target,
              pathway = pathway_id,
              pathway_name = pathway$name %||% pathway_id,
              pathway_type = pathway$type
            )
          }
        }
      }

      # Add output compound for this pathway (don't split, use full string)
      if (!is.null(pathway$output)) {
        compound_name <- pathway$output$compound
        kegg_id <- pathway$output$kegg_compound

        # Use KEGG ID for deduplication if available
        if (!is.null(kegg_id) && !is.na(kegg_id) && nchar(kegg_id) > 0) {
          compound_id <- paste0("COMPOUND_", kegg_id)
        } else {
          compound_id <- paste0("COMPOUND_", gsub("[^A-Za-z0-9]", "_", compound_name))
        }

        # Add compound info (don't override if already exists)
        if (!(compound_id %in% names(compound_info))) {
          compound_info[[compound_id]] <- list(
            id = compound_id,
            name = compound_name,
            kegg_id = kegg_id
          )
          compound_status[[compound_id]] <- "output"
        } else if (compound_status[[compound_id]] == "intermediate") {
          # Only upgrade from intermediate to output
          compound_status[[compound_id]] <- "output"
        }
        # If it's already input, keep it as input (green priority)

        # Connect source genes to output compound
        if (!is.null(pathway$output$sources)) {
          for (source in pathway$output$sources) {
            all_dag_edges[[length(all_dag_edges) + 1]] <- list(
              from = source,
              to = compound_id,
              pathway = pathway_id,
              pathway_name = pathway$name %||% pathway_id,
              pathway_type = pathway$type
            )
          }
        }
      }
    }
  } else {
    # Single-pathway potato (legacy schema)
    # Add input compound if present
    if (!is.null(potato@input)) {
      compound_info[["INPUT"]] <- list(
        id = "INPUT",
        name = potato@input$compound,
        from_step = 0,
        to_step = 1,
        kegg_id = potato@input$kegg_compound
      )
    }

    # Add output compound if present
    if (!is.null(potato@output)) {
      compound_info[["OUTPUT"]] <- list(
        id = "OUTPUT",
        name = potato@output$compound,
        from_step = 999,
        to_step = 1000,
        kegg_id = potato@output$kegg_compound
      )
    }

    for (edge in potato@edges) {
      # Skip edges with null endpoints (external compounds)
      if (is.null(edge$from) || is.null(edge$to)) {
        next
      }

      if (!is.null(edge$compound)) {
        # Extract steps from node IDs
        from_step <- as.integer(sub(".*_(\\d+)$", "\\1", edge$from))
        to_step <- as.integer(sub(".*_(\\d+)$", "\\1", edge$to))

        # Create compound node ID: compound_fromStep_toStep
        compound_id <- paste0("COMPOUND_", from_step, "_", to_step)
        compound_info[[compound_id]] <- list(
          id = compound_id,
          name = edge$compound,
          from_step = from_step,
          to_step = to_step,
          kegg_id = edge$kegg_compound
        )

        # Create edges: gene -> compound -> gene
        all_dag_edges[[length(all_dag_edges) + 1]] <- list(from = edge$from, to = compound_id)
        all_dag_edges[[length(all_dag_edges) + 1]] <- list(from = compound_id, to = edge$to)
      } else {
        # Direct gene -> gene edge
        all_dag_edges[[length(all_dag_edges) + 1]] <- list(from = edge$from, to = edge$to)
      }
    }

    # Add edges for INPUT/OUTPUT if present
    if (!is.null(potato@input) && !is.null(potato@input$targets)) {
      for (target in potato@input$targets) {
        all_dag_edges[[length(all_dag_edges) + 1]] <- list(from = "INPUT", to = target)
      }
    }

    if (!is.null(potato@output) && !is.null(potato@output$sources)) {
      for (source in potato@output$sources) {
        all_dag_edges[[length(all_dag_edges) + 1]] <- list(from = source, to = "OUTPUT")
      }
    }
  }

  # Build edge list for igraph
  if (length(all_dag_edges) == 0) {
    stop("No edges found in potato", call. = FALSE)
  }

  # Filter out edges with NULL, NA, or empty endpoints (external compounds)
  valid_edges <- Filter(function(e) {
    !is.null(e$from) && !is.null(e$to) &&
    length(e$from) > 0 && length(e$to) > 0 &&
    !is.na(e$from) && !is.na(e$to) &&
    nchar(as.character(e$from)) > 0 && nchar(as.character(e$to)) > 0
  }, all_dag_edges)

  if (length(valid_edges) == 0) {
    stop("No valid edges found in potato (all edges have NULL endpoints)", call. = FALSE)
  }

  edge_list <- do.call(rbind, lapply(valid_edges, function(e) {
    c(as.character(e$from), as.character(e$to))
  }))

  # Create graph
  g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)

  # Mark node types
  node_types <- sapply(igraph::V(g)$name, function(n) {
    if (grepl("^COMPOUND_", n)) {
      "compound"
    } else {
      "enzyme"
    }
  })
  igraph::V(g)$node_type <- node_types

  # Add compound names as vertex attribute
  compound_names <- sapply(igraph::V(g)$name, function(n) {
    if (n %in% names(compound_info)) {
      compound_info[[n]]$name
    } else {
      NA_character_
    }
  })
  igraph::V(g)$compound_name <- compound_names

  # Add compound KEGG IDs as vertex attribute
  compound_kegg_ids <- sapply(igraph::V(g)$name, function(n) {
    if (n %in% names(compound_info)) {
      kegg_id <- compound_info[[n]]$kegg_id
      if (!is.null(kegg_id) && !is.na(kegg_id)) kegg_id else NA_character_
    } else {
      NA_character_
    }
  })
  igraph::V(g)$kegg_id <- compound_kegg_ids

  # Add compound status (input, output, intermediate) as vertex attribute
  compound_status_vals <- sapply(igraph::V(g)$name, function(n) {
    if (n %in% names(compound_status)) {
      compound_status[[n]]
    } else {
      NA_character_
    }
  })
  igraph::V(g)$compound_status <- compound_status_vals

  g
}


#' Create step-based layout for single-pathway potatoes (internal)
#' @noRd
create_step_layout <- function(potato, node_names, is_bipartite = FALSE) {
  # Extract step numbers from node IDs (format: geneSymbol_stepNumber)
  steps <- sapply(node_names, function(n) {
    if (grepl("^COMPOUND_", n) || n %in% c("INPUT", "OUTPUT")) {
      # Handle compound nodes
      if (n == "INPUT") return(0)
      if (n == "OUTPUT") return(999)
      # Extract from_step from compound node ID
      parts <- strsplit(n, "_")[[1]]
      if (length(parts) >= 2) {
        return(as.integer(parts[2]))
      }
      return(NA)
    }
    # Extract step from gene node (last part after underscore)
    parts <- strsplit(n, "_")[[1]]
    as.integer(parts[length(parts)])
  })

  # Create layout: x = step number, y = staggered for alternatives at same step
  layout <- data.frame(name = node_names, x = steps, y = 0)

  # Stagger nodes at same step
  for (step in unique(steps[!is.na(steps)])) {
    step_nodes <- which(steps == step)
    if (length(step_nodes) > 1) {
      # Spread nodes vertically
      layout$y[step_nodes] <- seq(-1, 1, length.out = length(step_nodes))
    }
  }

  layout
}




#' Export potato DAG to graphviz DOT format
#'
#' @param potato Potato object
#' @param file Output file path (optional, returns string if NULL)
#' @export
export_potato_dot <- function(potato, file = NULL) {
  g <- build_potato_graph(potato)

  if (is.null(file)) {
    return(igraph::write_graph(g, file = tempfile(), format = "dot"))
  } else {
    igraph::write_graph(g, file = file, format = "dot")
  }
}
