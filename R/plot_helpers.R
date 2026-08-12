#' Prepare potato data for plotting (internal)
#' @noRd
prepare_potato_for_plotting <- function(potato, sack, genome_name, show_compounds, pathway) {

  # Load potato if character path
  if (is.character(potato)) {
    # Check schema version
    potato_data <- jsonlite::read_json(potato, simplifyVector = FALSE)
    if (!is.null(potato_data$schema_version) && potato_data$schema_version == "v2") {
      cli::cli_abort(c(
        "V2 schema potatoes are not yet supported by plot_potato_interactive2()",
        "i" = "Use plot_v2() instead:",
        " " = "pot <- load_potato_v2('{potato}')",
        " " = "g <- build_graph_v2(pot)",
        " " = "plot_v2(g, interactive = TRUE)"
      ))
    }
    potato <- load_potato(potato)
  }

  # Check if multi-pathway network
  is_multi_pathway <- !is.null(potato@edges) &&
                      is.list(potato@edges) &&
                      length(names(potato@edges)) > 0 &&
                      !is.null(potato@edges[[1]]$type)

  # Filter to pathway(s) if requested
  if (!is.null(pathway)) {
    if (!is_multi_pathway) {
      cli::cli_abort("The {.arg pathway} parameter only works with multi-pathway networks")
    }

    # Check all requested pathways exist
    missing <- setdiff(pathway, names(potato@edges))
    if (length(missing) > 0) {
      available <- paste(names(potato@edges), collapse = ", ")
      cli::cli_abort("Pathway(s) {.val {missing}} not found. Available: {available}")
    }

    potato_filtered <- potato
    potato_filtered@edges <- list()

    # Include all requested pathways
    all_node_ids <- character()
    for (pw in pathway) {
      potato_filtered@edges[[pw]] <- potato@edges[[pw]]
      all_node_ids <- c(all_node_ids, names(potato@edges[[pw]]$nodes))
    }

    # Get unique node IDs across all selected pathways
    all_node_ids <- unique(all_node_ids)
    potato_filtered@genes <- Filter(function(n) n$id %in% all_node_ids, potato@genes)

    potato <- potato_filtered
  }

  # Build graph
  if (show_compounds) {
    g <- build_bipartite_graph(potato)
  } else {
    g <- build_potato_graph(potato)
  }

  # Get node detection status
  # Check if genome_name provided without sack
  if (!is.null(genome_name) && is.null(sack)) {
    cli::cli_abort(c(
      "Cannot display genome {.val {genome_name}} without a sack",
      "i" = "Provide a {.arg sack} parameter with annotation results"
    ))
  }

  has_genome <- !is.null(sack) && !is.null(genome_name)
  if (has_genome) {
    # Validate genome exists in sack
    available_genomes <- sapply(sack@genomes, function(g) g$short_name)
    if (!genome_name %in% available_genomes) {
      cli::cli_abort(c(
        "Genome {.val {genome_name}} not found in sack",
        "i" = "Available genomes: {paste(available_genomes, collapse=', ')}"
      ))
    }
    node_status <- get_node_status(potato, sack, genome_name)
  } else {
    # Only include enzyme nodes initially - compound nodes will be added below
    enzyme_nodes <- igraph::V(g)$name[igraph::V(g)$node_type == "enzyme"]
    node_status <- tibble::tibble(
      name = enzyme_nodes,
      detected = NA,
      status = "Unknown"
    )
  }

  # Ensure all nodes in graph are in node_status
  all_nodes <- igraph::V(g)$name
  missing_nodes <- setdiff(all_nodes, node_status$name)

  if (length(missing_nodes) > 0) {
    # Get compound status from graph attributes (set during bipartite graph construction)
    compound_status_attr <- igraph::V(g)$compound_status

    missing_status <- tibble::tibble(
      name = missing_nodes,
      detected = NA,
      status = sapply(seq_along(missing_nodes), function(i) {
        n <- missing_nodes[i]
        # Use compound_status attribute if available
        idx <- which(igraph::V(g)$name == n)
        if (length(idx) > 0 && !is.null(compound_status_attr) && !is.na(compound_status_attr[idx])) {
          status <- compound_status_attr[idx]
          if (status == "input") return("Input")
          if (status == "output") return("Output")
          if (status == "intermediate") return("Compound")
        }
        # Fallback: infer from name prefix
        if (grepl("^COMPOUND_", n)) {
          "Compound"
        } else {
          "Unknown"
        }
      }),
      is_complex = FALSE,
      fraction_detected = NA
    )
    node_status <- dplyr::bind_rows(node_status, missing_status)
  }

  node_status <- node_status[match(all_nodes, node_status$name), ]

  # Add compound flag
  node_status$is_compound_node <- grepl("^COMPOUND_", node_status$name)

  # Gene IDs
  if (is_multi_pathway) {
    node_status$gene_id <- node_status$name
  } else {
    node_status$gene_id <- ifelse(node_status$is_compound_node,
                                   node_status$name,
                                   sub("_\\d+$", "", node_status$name))
  }

  # Get EC numbers, required, marker from potato genes
  node_status$ec <- purrr::map_chr(seq_len(nrow(node_status)), function(i) {
    if (node_status$is_compound_node[i]) return("")
    id <- node_status$gene_id[i]
    gene <- purrr::keep(potato@genes, ~ .x$id == id)
    if (length(gene) > 0 && !is.null(gene[[1]]$ec)) {
      ec_nums <- gene[[1]]$ec
      if (length(ec_nums) > 0) {
        return(paste0("\n[", ec_nums[1], "]"))
      }
    }
    return("")
  })

  # Get required and marker status (pathway-specific in multi-pathway networks)
  node_status$required <- purrr::map_lgl(seq_len(nrow(node_status)), function(i) {
    if (node_status$is_compound_node[i]) return(FALSE)
    id <- node_status$gene_id[i]

    # For multi-pathway networks, check pathways
    if (is_multi_pathway) {
      # Gene is required if required in ANY pathway
      for (pathway_id in names(potato@edges)) {
        pathway <- potato@edges[[pathway_id]]
        if (id %in% names(pathway$nodes)) {
          pathway_node <- pathway$nodes[[id]]
          if (!is.null(pathway_node$required) && pathway_node$required) {
            return(TRUE)
          }
        }
      }
      return(FALSE)
    } else {
      # Single-pathway potato
      gene <- purrr::keep(potato@genes, ~ .x$id == id)
      if (length(gene) > 0) return(gene[[1]]$required %||% FALSE)
      return(FALSE)
    }
  })

  node_status$marker <- purrr::map_lgl(seq_len(nrow(node_status)), function(i) {
    if (node_status$is_compound_node[i]) return(FALSE)
    id <- node_status$gene_id[i]

    # For multi-pathway networks, check pathways
    if (is_multi_pathway) {
      # Gene is marker if marker in ANY pathway
      for (pathway_id in names(potato@edges)) {
        pathway <- potato@edges[[pathway_id]]
        if (id %in% names(pathway$nodes)) {
          pathway_node <- pathway$nodes[[id]]
          if (!is.null(pathway_node$marker) && pathway_node$marker) {
            return(TRUE)
          }
        }
      }
      return(FALSE)
    } else {
      # Single-pathway potato
      gene <- purrr::keep(potato@genes, ~ .x$id == id)
      if (length(gene) > 0) return(gene[[1]]$marker %||% FALSE)
      return(FALSE)
    }
  })

  # Find pathway membership for multi-pathway networks
  node_pathways <- character(nrow(node_status))
  if (is_multi_pathway) {
    # Use graph attribute if available
    if (!is.null(igraph::V(g)$pathway_membership)) {
      node_pathways <- igraph::V(g)$pathway_membership[match(node_status$name, igraph::V(g)$name)]
      node_pathways[is.na(node_pathways)] <- ""
    } else {
      # Look up by gene ID
      for (i in seq_len(nrow(node_status))) {
        node_id <- node_status$gene_id[i]
        pathways_with_node <- character()
        for (pathway_id in names(potato@edges)) {
          pathway <- potato@edges[[pathway_id]]
          if (node_id %in% names(pathway$nodes)) {
            pathways_with_node <- c(pathways_with_node, pathway$name %||% pathway_id)
          }
        }
        node_pathways[i] <- paste(pathways_with_node, collapse = ", ")
      }
    }
  }

  # Build hover text (HTML format for visNetwork)
  hover_texts <- character(nrow(node_status))
  for (i in seq_len(nrow(node_status))) {
    if (node_status$is_compound_node[i]) {
      cmp_name <- igraph::V(g)$compound_name[match(node_status$name[i], igraph::V(g)$name)]
      kegg_id <- igraph::V(g)$kegg_id[match(node_status$name[i], igraph::V(g)$name)]
      hover_texts[i] <- paste0("Compound: ", cmp_name)
      if (!is.null(kegg_id) && !is.na(kegg_id) && nchar(kegg_id) > 0) {
        hover_texts[i] <- paste0(hover_texts[i], "<br>KEGG: ", kegg_id)
      }
    } else {
      ec_clean <- gsub("^\n\\[|\\]$", "", node_status$ec[i])
      type_str <- ifelse(node_status$marker[i], "Marker", ifelse(node_status$required[i], "Required", "Optional"))

      # Get databases info
      id <- node_status$gene_id[i]
      gene <- purrr::keep(potato@genes, ~ .x$id == id)
      databases_str <- ""
      if (length(gene) > 0 && !is.null(gene[[1]]$databases)) {
        dbs <- gene[[1]]$databases
        db_lines <- character()
        for (db_name in names(dbs)) {
          terms <- paste(dbs[[db_name]], collapse = ", ")
          db_lines <- c(db_lines, paste0(db_name, ": ", terms))
        }
        if (length(db_lines) > 0) {
          databases_str <- paste0("<br>", paste(db_lines, collapse = "<br>"))
        }
      }

      # Format pathway list with newlines
      pathway_str <- ""
      if (is_multi_pathway && node_pathways[i] != "") {
        pathway_list <- strsplit(node_pathways[i], ", ")[[1]]
        pathway_lines <- paste0(" - ", paste(pathway_list, collapse = "<br> - "))
        pathway_str <- paste0("Pathways:<br>", pathway_lines, "<br>")
      }

      hover_texts[i] <- paste0(
        "Gene: ", node_status$gene_id[i],
        if (node_status$ec[i] != "") paste0(" EC: ", ec_clean) else "", "<br>",
        "Type: ", type_str, "<br>",
        pathway_str,
        "Status: ", node_status$status[i],
        databases_str
      )
    }
  }

  # Add shape and label
  node_status <- node_status %>%
    dplyr::mutate(
      node_shape = dplyr::case_when(
        is_compound_node ~ "triangle",
        TRUE ~ "circle"  # All genes are circles now (no diamond for markers)
      ),
      label = ifelse(
        is_compound_node,
        igraph::V(g)$compound_name[match(name, igraph::V(g)$name)],
        paste0(gene_id, ec)
      ),
      hover_text = hover_texts
    )

  list(
    potato = potato,
    g = g,
    node_status = node_status,
    is_multi_pathway = is_multi_pathway,
    has_genome = has_genome
  )
}


#' Calculate node layout coordinates (internal)
#' @noRd
calculate_node_layout <- function(potato, g, is_multi_pathway, show_compounds, layout = "xy", scale_for_visnetwork = FALSE) {

  has_curated_coords <- FALSE
  use_curated <- (layout == "curated" || layout == "xy")  # Use curated coords for "curated" or legacy "xy"

  # Dynamic scaling based on network size
  # visNetwork needs larger coordinates for bigger networks
  # Base scale increases with number of nodes (sqrt to avoid extreme scaling)
  if (scale_for_visnetwork) {
    n_nodes <- igraph::vcount(g)
    # Scale from ~20 (small networks) to ~60 (large networks)
    coord_scale <- 20 + (sqrt(n_nodes) * 2)
    coord_scale <- max(20, min(coord_scale, 60))  # Clamp between 20-60
  } else {
    coord_scale <- 1
  }

  if (is_multi_pathway && use_curated) {
    # Check for curated coordinates
    x_field <- if (show_compounds) "x_compounds" else "x"
    y_field <- if (show_compounds) "y_compounds" else "y"

    node_coords_list <- list()
    for (gene in potato@genes) {
      x_val <- gene[[x_field]]
      y_val <- gene[[y_field]]

      # Fall back to regular x,y if compound-specific coords don't exist
      if (show_compounds && (is.null(x_val) || is.null(y_val))) {
        x_val <- gene$x
        y_val <- gene$y
      }

      if (!is.null(x_val) && !is.null(y_val)) {
        node_coords_list[[gene$id]] <- c(x_val, y_val)
      }
    }

    # If we have coordinates for some enzyme nodes, use hybrid approach
    enzyme_nodes <- igraph::V(g)$name[igraph::V(g)$node_type == "enzyme"]
    nodes_with_coords <- enzyme_nodes[enzyme_nodes %in% names(node_coords_list)]

    if (length(nodes_with_coords) > 0) {
      has_curated_coords <- TRUE

      # Start with layout algorithm for all nodes
      # For "xy", use fr as base for nodes without curated coords
      layout_matrix <- switch(layout,
        "xy" = if (show_compounds) igraph::layout_with_kk(g) else igraph::layout_with_fr(g),
        "sugiyama" = igraph::layout_with_sugiyama(g)$layout,
        "fr" = igraph::layout_with_fr(g),
        "kk" = igraph::layout_with_kk(g),
        "tree" = igraph::layout_as_tree(g),
        "circle" = igraph::layout_in_circle(g),
        "grid" = igraph::layout_on_grid(g),
        if (show_compounds) igraph::layout_with_kk(g) else igraph::layout_with_fr(g)
      )

      node_coords <- data.frame(
        name = igraph::V(g)$name,
        x = layout_matrix[, 1],
        y = layout_matrix[, 2]
      )

      # Override with curated coordinates where available
      for (i in seq_along(node_coords$name)) {
        node_name <- node_coords$name[i]
        if (node_name %in% names(node_coords_list)) {
          node_coords$x[i] <- node_coords_list[[node_name]][1] * coord_scale
          node_coords$y[i] <- node_coords_list[[node_name]][2] * coord_scale
        }
      }

      # For compound nodes, use stored coordinates or position between connected enzymes
      if (show_compounds) {
        # Check if potato has compound_coordinates
        compound_coords_stored <- potato@compound_coordinates %||% list()

        compound_nodes <- which(igraph::V(g)$node_type == "compound")
        for (i in compound_nodes) {
          compound_id <- node_coords$name[i]

          # Use stored coordinates if available
          if (compound_id %in% names(compound_coords_stored)) {
            stored <- compound_coords_stored[[compound_id]]
            node_coords$x[i] <- stored$x * coord_scale
            node_coords$y[i] <- stored$y * coord_scale
          } else {
            # Calculate from neighbors
            neighbors <- igraph::neighbors(g, i, mode = "all")
            neighbor_names <- igraph::V(g)$name[neighbors]
            neighbor_coords <- node_coords[node_coords$name %in% neighbor_names &
                                           node_coords$name %in% names(node_coords_list), ]
            if (nrow(neighbor_coords) > 0) {
              node_coords$x[i] <- mean(neighbor_coords$x)
              node_coords$y[i] <- mean(neighbor_coords$y)
            }
          }
        }
      }

      cli::cli_alert_info("Using curated coordinates for {length(nodes_with_coords)}/{length(enzyme_nodes)} nodes")
    }
  }

  # If no curated coordinates, use layout algorithms
  if (!has_curated_coords) {
    if (is_multi_pathway) {
      # If layout="xy" but no curated coords, fall back to "fr"
      actual_layout <- if (layout == "xy") "fr" else layout

      layout_matrix <- switch(actual_layout,
        "sugiyama" = igraph::layout_with_sugiyama(g)$layout,
        "fr" = igraph::layout_with_fr(g),
        "kk" = igraph::layout_with_kk(g),
        "tree" = igraph::layout_as_tree(g),
        "circle" = igraph::layout_in_circle(g),
        "grid" = igraph::layout_on_grid(g),
        if (show_compounds) igraph::layout_with_kk(g) else igraph::layout_with_fr(g)
      )
      node_coords <- data.frame(
        name = igraph::V(g)$name,
        x = layout_matrix[, 1],
        y = layout_matrix[, 2]
      )
    } else {
      # Single-pathway uses step-based layout
      node_coords <- create_step_layout(potato, igraph::V(g)$name, is_bipartite = show_compounds)
    }
  }

  node_coords
}
