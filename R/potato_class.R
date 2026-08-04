#' Potato S7 Class
#'
#' Represents a pathway definition with genes and pathway structure
#'
#' @param id Character. Unique identifier for the potato
#' @param name Character. Human-readable name
#' @param nodes List. Nodes (genes/enzymes) in the pathway. Each node can have:
#'   - `marker`: Logical. If TRUE, this gene is a diagnostic marker for the pathway
#'   - `required`: Logical. If TRUE, pathway cannot be complete without this gene
#'   - `type`: Character. "enzyme", "compound", or "transporter"
#'   - `databases`: List. Detection methods by database type (e.g., kofam, blast, hmm)
#' @param edges List. Edges connecting nodes
#' @param tags Character vector. Tags for organizing potatoes
#' @param source Character. Source/origin of the pathway definition
#' @param notes Character. Additional notes about the pathway
#' @param scoring List. Scoring parameters including:
#'   - `min_fraction`: Minimum fraction of genes required (default 0.75)
#'   - `marker_mode`: "any" or "all" - how many marker genes needed for positive call
#' @param json_path Character. Path to the source JSON file
#' @param graph igraph object or NULL. Cached graph representation
#'
#' @export
Potato <- S7::new_class(
  "Potato",
  properties = list(
    id = S7::class_character,
    name = S7::class_character,
    nodes = S7::class_list,
    edges = S7::class_list,
    tags = S7::class_character,
    source = S7::class_character,
    notes = S7::class_character,
    scoring = S7::class_list,
    input = S7::class_list,
    output = S7::class_list,
    json_path = S7::class_character,
    graph = S7::class_any,
    compound_coordinates = S7::class_list
  ),
  validator = function(self) {
    # Basic validation
    if (nchar(self@id) == 0) {
      "Potato must have an id"
    } else if (nchar(self@name) == 0) {
      "Potato must have a name"
    } else if (length(self@nodes) == 0) {
      "Potato must have at least one node"
    }
  }
)


#' Load a potato from JSON
#'
#' Loads multi-pathway network potatoes (new schema). Does not support
#' single-pathway schema.
#'
#' @param path Path to potato JSON file
#' @return Potato object (for multi-pathway networks, stores pathways in edges slot)
#' @export
load_potato <- function(path) {
  if (!file.exists(path)) {
    stop("Potato file not found: ", path, call. = FALSE)
  }

  data <- jsonlite::read_json(path, simplifyVector = FALSE)

  # Check if single-pathway potato
  is_network <- !is.null(data$pathways) && is.list(data$pathways)

  if (!is_network) {
    # Check if marked as inactive
    if (!is.null(data$active) && data$active == FALSE) {
      warning("Potato '", data$id, "' is inactive (active: false). ",
              "Consider using the updated multi-pathway network version.",
              call. = FALSE)
    } else {
      stop("Single-pathway schema is no longer supported. ",
           "Only multi-pathway network potatoes can be loaded.",
           call. = FALSE)
    }
  }

  # Load multi-pathway network
  nodes <- if (is.null(data$nodes)) list() else data$nodes
  tags <- if (is.null(data$tags)) character(0) else unlist(data$tags)

  # Store pathways in edges slot (temporary until Potato class updated)
  # This maintains compatibility with existing S7 class structure
  edges <- if (is.null(data$pathways)) list() else data$pathways

  Potato(
    id = data$id,
    name = data$name,
    nodes = nodes,
    edges = edges,  # Contains pathways for network potatoes
    tags = tags,
    source = if (is.null(data$source)) "" else data$source,
    notes = if (is.null(data$notes)) "" else data$notes,
    scoring = if (is.null(data$scoring)) list() else data$scoring,
    input = if (is.null(data$input)) list() else data$input,
    output = if (is.null(data$output)) list() else data$output,
    json_path = path,
    graph = NULL,
    compound_coordinates = if (is.null(data$compound_coordinates)) list() else data$compound_coordinates
  )
}


#' Load multiple potatoes from directory
#'
#' By default, only loads active potatoes (active != false). Set include_inactive = TRUE
#' to load inactive potatoes.
#'
#' @param dir Directory containing potato JSON files
#' @param tags Optional character vector of tags to filter by
#' @param include_inactive Logical. If TRUE, loads inactive potatoes (active: false)
#' @return Named list of Potato objects
#' @export
load_potatoes <- function(dir, tags = NULL, include_inactive = FALSE) {
  if (!dir.exists(dir)) {
    stop("Directory not found: ", dir, call. = FALSE)
  }

  json_files <- list.files(dir, pattern = "\\.json$", full.names = TRUE)

  potatoes <- lapply(json_files, function(filepath) {
    tryCatch({
      # Read JSON to check active flag before loading
      data <- jsonlite::read_json(filepath, simplifyVector = FALSE)

      # Skip inactive potatoes unless requested
      if (!include_inactive && !is.null(data$active) && data$active == FALSE) {
        return(NULL)
      }

      # Skip single-pathway potatoes that aren't marked inactive
      is_network <- !is.null(data$pathways) && is.list(data$pathways)
      if (!is_network && (is.null(data$active) || data$active == TRUE)) {
        warning("Skipping ", basename(filepath), ": single-pathway schema not supported (should be marked active: false)", call. = FALSE)
        return(NULL)
      }

      # Load the potato
      load_potato(filepath)
    }, error = function(e) {
      warning("Failed to load ", basename(filepath), ": ", e$message, call. = FALSE)
      NULL
    })
  })

  # Remove NULLs from skipped/failed loads
  potatoes <- potatoes[!sapply(potatoes, is.null)]

  if (length(potatoes) == 0) {
    message("No active potatoes found in ", dir)
    return(list())
  }

  # Name by ID
  names(potatoes) <- sapply(potatoes, function(p) p@id)

  # Filter by tags if specified
  if (!is.null(tags)) {
    potatoes <- potatoes[sapply(potatoes, function(p) {
      any(tags %in% p@tags)
    })]
  }

  potatoes
}


#' Load example test potato
#'
#' @return Potato object
#' @export
load_test_potato <- function() {
  test_path <- system.file("potatoes", "test_glycolysis.json", package = "potato")
  if (test_path == "") {
    stop("Test potato not found in package", call. = FALSE)
  }
  load_potato(test_path)
}


#' Get enzyme nodes from potato
#'
#' @param potato Potato object
#' @return List of enzyme nodes
#' @keywords internal
get_enzyme_nodes <- function(potato) {
  S7::check_is_S7(potato)

  # Filter nodes that are enzymes
  enzyme_nodes <- Filter(function(node) {
    !is.null(node$type) && node$type == "enzyme"
  }, potato@nodes)

  enzyme_nodes
}


#' Get detection terms for a database type
#'
#' Extracts all terms (KO IDs, gene names, etc.) for a specific database type from
#' a potato's enzyme nodes.
#'
#' @param potato Potato object
#' @param database_name Database type ("kofam", "blast", "hmm", "pfam")
#' @return Character vector of unique terms
#' @keywords internal
get_detection_terms <- function(potato, database_name) {
  S7::check_is_S7(potato)

  enzyme_nodes <- get_enzyme_nodes(potato)

  # Extract terms from databases field
  terms <- unlist(lapply(enzyme_nodes, function(node) {
    if (!is.null(node$databases) && !is.null(node$databases[[database_name]])) {
      node$databases[[database_name]]
    }
  }))

  unique(terms)
}


#' Get marker gene nodes from potato
#'
#' @param potato Potato object
#' @return List of nodes marked as marker genes
#' @keywords internal
get_marker_genes <- function(potato) {
  S7::check_is_S7(potato)

  # Filter nodes that are marked as markers
  marker_nodes <- Filter(function(node) {
    !is.null(node$marker) && node$marker == TRUE
  }, potato@nodes)

  marker_nodes
}


#' Build igraph from potato
#'
#' @param potato Potato object
#' @return igraph object
#' @keywords internal
build_potato_graph <- function(potato) {
  S7::check_is_S7(potato)

  # Check if this is a multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  if (is_network) {
    # Multi-pathway network: create gene-based graph (not step-based)
    # Genes are shared across pathways - same gene appears once
    edge_pathways <- list()  # Key: "from_gene|to_gene", Value: list of pathway info

    for (pathway_id in names(potato@edges)) {
      pathway <- potato@edges[[pathway_id]]

      if (is.null(pathway$edges) || length(pathway$edges) == 0) {
        next
      }

      # For multi-pathway networks, use gene IDs directly (not id_step)
      # This creates a connected graph where genes are shared
      for (edge in pathway$edges) {
        from_id <- edge$from
        to_id <- edge$to

        edge_key <- paste0(from_id, "|", to_id)

        # Add this pathway to the edge's pathway list
        if (is.null(edge_pathways[[edge_key]])) {
          edge_pathways[[edge_key]] <- list(
            from = from_id,
            to = to_id,
            pathways = character(),
            pathway_names = character(),
            pathway_types = character()
          )
        }

        edge_pathways[[edge_key]]$pathways <- c(edge_pathways[[edge_key]]$pathways, pathway_id)
        edge_pathways[[edge_key]]$pathway_names <- c(edge_pathways[[edge_key]]$pathway_names, pathway$name %||% pathway_id)
        edge_pathways[[edge_key]]$pathway_types <- c(edge_pathways[[edge_key]]$pathway_types, pathway$type)
      }
    }

    if (length(edge_pathways) == 0) {
      # No edges in any pathway - create graph with just nodes
      # Collect all unique gene IDs
      all_gene_ids <- character()
      for (pathway_id in names(potato@edges)) {
        pathway <- potato@edges[[pathway_id]]
        all_gene_ids <- c(all_gene_ids, names(pathway$nodes))
      }
      all_gene_ids <- unique(all_gene_ids)

      g <- igraph::make_empty_graph(n = length(all_gene_ids), directed = TRUE)
      igraph::V(g)$name <- all_gene_ids
      igraph::V(g)$node_type <- "enzyme"
      return(g)
    }

    # Build edge list from unique edges
    edge_list <- do.call(rbind, lapply(edge_pathways, function(e) {
      c(e$from, e$to)
    }))

    g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)

    # Add edge attributes for pathway membership
    igraph::E(g)$pathways <- sapply(edge_pathways, function(e) paste(e$pathways, collapse = ","))
    igraph::E(g)$pathway_names <- sapply(edge_pathways, function(e) paste(e$pathway_names, collapse = ","))
    igraph::E(g)$n_pathways <- sapply(edge_pathways, function(e) length(e$pathways))
    igraph::E(g)$is_shared <- sapply(edge_pathways, function(e) length(e$pathways) > 1)

    # Add node attributes for pathway membership and type
    node_pathway_membership <- sapply(igraph::V(g)$name, function(gene_id) {
      pathways <- character()
      for (pathway_id in names(potato@edges)) {
        if (gene_id %in% names(potato@edges[[pathway_id]]$nodes)) {
          pathways <- c(pathways, potato@edges[[pathway_id]]$name %||% pathway_id)
        }
      }
      paste(pathways, collapse = ", ")
    })
    igraph::V(g)$pathway_membership <- node_pathway_membership

    # Mark all nodes as enzymes (multi-pathway networks don't have compound nodes in the graph)
    igraph::V(g)$node_type <- "enzyme"

    return(g)

  } else {
    # Single-pathway potato (legacy schema)
    if (length(potato@edges) == 0) {
      # Graph with just nodes, no edges
      # Use the 'nodes' field which has id_step format, not just 'id'
      node_names <- unlist(sapply(potato@nodes, function(n) n$nodes))
      g <- igraph::make_empty_graph(n = length(node_names), directed = TRUE)
      igraph::V(g)$name <- node_names
      return(g)
    }

    # Build edge list
    edge_list <- do.call(rbind, lapply(potato@edges, function(e) {
      c(e$from, e$to)
    }))

    # Create graph
    g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)
    return(g)
  }
}


#' Validate potato structure
#'
#' Validates a potato JSON structure for common errors and required fields.
#' Handles both multi-pathway networks (new schema) and single-pathway potatoes (legacy schema).
#' Users should run this on custom potatoes before using them for annotation.
#'
#' @param potato Potato object or list (raw JSON data)
#' @param strict Logical. If TRUE, performs additional strict checks (default FALSE)
#' @return List with 'valid' (logical), 'errors' (character vector), and 'warnings' (character vector)
#' @export
validate_potato <- function(potato, strict = FALSE) {
  # Handle both Potato S7 objects and raw JSON lists
  # S7 objects have class c("potato::Potato", "S7_object")
  is_s7_potato <- "S7_object" %in% class(potato)
  is_list_with_id <- is.list(potato) && "id" %in% names(potato)

  if (is_s7_potato || is_list_with_id) {
    if (is_s7_potato) {
      # Extract data from S7 object
      # NOTE: For multi-pathway networks, edges slot contains pathways
      data <- list(
        id = potato@id,
        name = potato@name,
        nodes = potato@nodes,
        tags = potato@tags,
        source = potato@source,
        scoring = potato@scoring
      )

      # Check if edges slot contains pathways (multi-pathway network)
      # Pathways have names and first element has a 'type' field
      if (!is.null(potato@edges) &&
          is.list(potato@edges) &&
          length(names(potato@edges)) > 0 &&
          !is.null(potato@edges[[1]]$type)) {
        data$pathways <- potato@edges
      } else {
        data$edges <- potato@edges
      }
    } else {
      data <- potato
    }
  } else {
    stop("potato must be a Potato object or list", call. = FALSE)
  }

  errors <- character(0)
  warnings <- character(0)

  # Check required fields
  if (is.null(data$id) || nchar(data$id) == 0) {
    errors <- c(errors, "Missing required field: 'id'")
  }
  if (is.null(data$name) || nchar(data$name) == 0) {
    errors <- c(errors, "Missing required field: 'name'")
  }

  # Check ID format
  if (!is.null(data$id) && !grepl("^[a-zA-Z0-9_-]+$", data$id)) {
    warnings <- c(warnings, "Potato ID should contain only letters, numbers, underscores, and hyphens")
  }

  # Detect if this is a multi-pathway network
  is_network <- !is.null(data$pathways) && is.list(data$pathways)

  if (is_network) {
    # Validate multi-pathway network schema
    result <- validate_multi_pathway(data, strict)
    return(result)
  } else {
    # Validate single-pathway schema (legacy format)
    result <- validate_single_pathway(data, strict)
    return(result)
  }
}

#' Validate single-pathway potato (legacy schema)
#' @keywords internal
validate_single_pathway <- function(data, strict) {
  errors <- character(0)
  warnings <- character(0)

  # Check nodes
  if (is.null(data$nodes) || length(data$nodes) == 0) {
    errors <- c(errors, "Potato has no nodes")
  }

  # Check node IDs are unique
  node_ids <- sapply(data$nodes, function(n) n$id)
  if (length(node_ids) != length(unique(node_ids))) {
    errors <- c(errors, "Duplicate node IDs found")
  }

  # Validate each node
  for (i in seq_along(data$nodes)) {
    node <- data$nodes[[i]]
    node_prefix <- sprintf("Node %d (%s)", i, node$id %||% "unnamed")

    # Required node fields
    if (is.null(node$id)) {
      errors <- c(errors, sprintf("%s: missing 'id'", node_prefix))
    }
    if (is.null(node$type)) {
      errors <- c(errors, sprintf("%s: missing 'type'", node_prefix))
    } else if (!node$type %in% c("enzyme", "compound", "transporter")) {
      warnings <- c(warnings, sprintf("%s: type '%s' is non-standard", node_prefix, node$type))
    }

    # Check enzyme nodes have detection methods
    if (!is.null(node$type) && node$type == "enzyme") {
      has_databases <- !is.null(node$databases) && length(node$databases) > 0

      if (!has_databases) {
        errors <- c(errors, sprintf("%s: enzyme missing 'databases' field", node_prefix))
      }

      # Validate databases field
      if (has_databases) {
        if (!is.list(node$databases)) {
          errors <- c(errors, sprintf("%s: 'databases' must be a list", node_prefix))
        } else {
          valid_db_types <- c("kofam", "blast", "hmm", "patric")
          for (db_name in names(node$databases)) {
            if (!db_name %in% valid_db_types) {
              errors <- c(errors, sprintf("%s: invalid database type '%s'", node_prefix, db_name))
            }

            db_terms <- node$databases[[db_name]]
            if (!is.character(db_terms) && !is.list(db_terms)) {
              errors <- c(errors, sprintf("%s: database '%s' terms must be character vector", node_prefix, db_name))
            }

            # Validate KO format
            if (db_name == "kofam") {
              ko_ids <- if (is.list(db_terms)) unlist(db_terms) else db_terms
              invalid_ko <- ko_ids[!grepl("^K[0-9]{5}$", ko_ids)]
              if (length(invalid_ko) > 0) {
                warnings <- c(warnings, sprintf("%s: invalid KO format: %s", node_prefix, paste(invalid_ko, collapse = ", ")))
              }
            }
          }
        }
      }

      # Check field types
      if (!is.null(node$required) && !is.logical(node$required)) {
        errors <- c(errors, sprintf("%s: 'required' must be TRUE/FALSE", node_prefix))
      }
      if (!is.null(node$marker) && !is.logical(node$marker)) {
        errors <- c(errors, sprintf("%s: 'marker' must be TRUE/FALSE", node_prefix))
      }
    }
  }

  # Collect valid node IDs
  all_valid_ids <- node_ids
  for (node in data$nodes) {
    if (!is.null(node$nodes)) {
      dag_node_ids <- if (is.list(node$nodes)) unlist(node$nodes) else node$nodes
      all_valid_ids <- c(all_valid_ids, dag_node_ids)
    }
  }
  all_valid_ids <- unique(all_valid_ids)

  # Validate edges
  if (!is.null(data$edges)) {
    for (i in seq_along(data$edges)) {
      edge <- data$edges[[i]]
      if (is.null(edge$from)) {
        errors <- c(errors, sprintf("Edge %d: missing 'from'", i))
      } else if (!edge$from %in% all_valid_ids) {
        errors <- c(errors, sprintf("Edge %d: 'from' references non-existent node '%s'", i, edge$from))
      }
      if (is.null(edge$to)) {
        errors <- c(errors, sprintf("Edge %d: missing 'to'", i))
      } else if (!edge$to %in% all_valid_ids) {
        errors <- c(errors, sprintf("Edge %d: 'to' references non-existent node '%s'", i, edge$to))
      }
    }
  }

  # Check for cycles
  if (length(errors) == 0 && !is.null(data$edges) && length(data$edges) > 0) {
    tryCatch({
      edge_list <- do.call(rbind, lapply(data$edges, function(e) c(e$from, e$to)))
      g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)
      if (!igraph::is_dag(g)) {
        errors <- c(errors, "Potato contains cycles (must be DAG)")
      }
    }, error = function(e) {
      errors <- c(errors, sprintf("Graph validation failed: %s", e$message))
    })
  }

  # Check marker genes
  if (!is.null(data$scoring$marker_mode)) {
    marker_nodes <- Filter(function(n) !is.null(n$marker) && n$marker == TRUE, data$nodes)
    if (length(marker_nodes) == 0) {
      warnings <- c(warnings, "marker_mode specified but no marker genes found")
    }
    if (!data$scoring$marker_mode %in% c("any", "all")) {
      errors <- c(errors, "marker_mode must be 'any' or 'all'")
    }
  }

  # Strict checks
  if (strict) {
    if (is.null(data$source) || nchar(data$source) == 0) {
      warnings <- c(warnings, "No source specified")
    }
    if (is.null(data$tags) || length(data$tags) == 0) {
      warnings <- c(warnings, "No tags specified")
    }
  }

  list(valid = length(errors) == 0, errors = errors, warnings = warnings)
}

#' Validate multi-pathway network potato (new schema)
#' @keywords internal
validate_multi_pathway <- function(data, strict) {
  errors <- character(0)
  warnings <- character(0)

  # Check nodes array exists
  if (is.null(data$nodes) || length(data$nodes) == 0) {
    errors <- c(errors, "Multi-pathway network has no global nodes")
  }

  # Check node IDs are unique
  node_ids <- sapply(data$nodes, function(n) n$id)
  if (length(node_ids) != length(unique(node_ids))) {
    errors <- c(errors, "Duplicate node IDs in global nodes array")
  }

  # Validate global nodes (detection methods only, no step/required/marker)
  for (i in seq_along(data$nodes)) {
    node <- data$nodes[[i]]
    node_prefix <- sprintf("Global node %d (%s)", i, node$id %||% "unnamed")

    if (is.null(node$id)) {
      errors <- c(errors, sprintf("%s: missing 'id'", node_prefix))
    }
    if (is.null(node$type)) {
      errors <- c(errors, sprintf("%s: missing 'type'", node_prefix))
    }

    # Global nodes should NOT have pathway-specific attributes
    if (!is.null(node$step)) {
      errors <- c(errors, sprintf("%s: 'step' belongs in pathway-specific nodes, not global", node_prefix))
    }
    if (!is.null(node$required)) {
      errors <- c(errors, sprintf("%s: 'required' belongs in pathway-specific nodes, not global", node_prefix))
    }
    if (!is.null(node$marker)) {
      errors <- c(errors, sprintf("%s: 'marker' belongs in pathway-specific nodes, not global", node_prefix))
    }

    # Validate databases field
    if (!is.null(node$type) && node$type == "enzyme") {
      if (is.null(node$databases) || length(node$databases) == 0) {
        errors <- c(errors, sprintf("%s: enzyme missing 'databases' field", node_prefix))
      } else {
        valid_db_types <- c("kofam", "blast", "hmm", "patric")
        for (db_name in names(node$databases)) {
          if (!db_name %in% valid_db_types) {
            errors <- c(errors, sprintf("%s: invalid database type '%s'", node_prefix, db_name))
          }

          # Validate KO format
          if (db_name == "kofam") {
            ko_ids <- unlist(node$databases[[db_name]])
            invalid_ko <- ko_ids[!grepl("^K[0-9]{5}$", ko_ids)]
            if (length(invalid_ko) > 0) {
              warnings <- c(warnings, sprintf("%s: invalid KO format: %s", node_prefix, paste(invalid_ko, collapse = ", ")))
            }
          }
        }
      }
    }
  }

  # Validate pathways field
  if (is.null(data$pathways) || !is.list(data$pathways) || length(data$pathways) == 0) {
    errors <- c(errors, "Multi-pathway network must have 'pathways' field with at least one pathway")
  } else {
    # Validate each pathway
    for (pathway_id in names(data$pathways)) {
      pathway <- data$pathways[[pathway_id]]
      path_prefix <- sprintf("Pathway '%s'", pathway_id)

      # Check pathway has nodes
      if (is.null(pathway$nodes) || length(pathway$nodes) == 0) {
        errors <- c(errors, sprintf("%s: no nodes defined", path_prefix))
      }

      # Check pathway type
      if (is.null(pathway$type)) {
        warnings <- c(warnings, sprintf("%s: missing 'type' field (should be 'variant' or 'independent')", path_prefix))
      } else if (!pathway$type %in% c("variant", "independent")) {
        errors <- c(errors, sprintf("%s: type must be 'variant' or 'independent', got '%s'", path_prefix, pathway$type))
      }

      # Validate pathway-specific node references
      if (!is.null(pathway$nodes)) {
        for (node_id in names(pathway$nodes)) {
          if (!node_id %in% node_ids) {
            errors <- c(errors, sprintf("%s: node '%s' not found in global nodes array", path_prefix, node_id))
          }

          node_attrs <- pathway$nodes[[node_id]]

          # Check required attributes
          if (is.null(node_attrs$step)) {
            errors <- c(errors, sprintf("%s node '%s': missing 'step'", path_prefix, node_id))
          }
          if (is.null(node_attrs$required)) {
            warnings <- c(warnings, sprintf("%s node '%s': missing 'required' (recommended)", path_prefix, node_id))
          }
          if (is.null(node_attrs$marker)) {
            warnings <- c(warnings, sprintf("%s node '%s': missing 'marker' (recommended)", path_prefix, node_id))
          }

          # Validate types
          if (!is.null(node_attrs$required) && !is.logical(node_attrs$required)) {
            errors <- c(errors, sprintf("%s node '%s': 'required' must be TRUE/FALSE", path_prefix, node_id))
          }
          if (!is.null(node_attrs$marker) && !is.logical(node_attrs$marker)) {
            errors <- c(errors, sprintf("%s node '%s': 'marker' must be TRUE/FALSE", path_prefix, node_id))
          }
        }
      }

      # Validate pathway edges
      if (!is.null(pathway$edges)) {
        pathway_node_ids <- names(pathway$nodes)
        for (i in seq_along(pathway$edges)) {
          edge <- pathway$edges[[i]]
          if (is.null(edge$from)) {
            errors <- c(errors, sprintf("%s edge %d: missing 'from'", path_prefix, i))
          } else if (!edge$from %in% pathway_node_ids) {
            errors <- c(errors, sprintf("%s edge %d: 'from' node '%s' not in pathway nodes", path_prefix, i, edge$from))
          }
          if (is.null(edge$to)) {
            errors <- c(errors, sprintf("%s edge %d: missing 'to'", path_prefix, i))
          } else if (!edge$to %in% pathway_node_ids) {
            errors <- c(errors, sprintf("%s edge %d: 'to' node '%s' not in pathway nodes", path_prefix, i, edge$to))
          }
        }

        # Check for cycles in this pathway
        if (length(errors) == 0) {
          tryCatch({
            edge_list <- do.call(rbind, lapply(pathway$edges, function(e) c(e$from, e$to)))
            g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)
            if (!igraph::is_dag(g)) {
              errors <- c(errors, sprintf("%s: contains cycles (must be DAG)", path_prefix))
            }
          }, error = function(e) {
            errors <- c(errors, sprintf("%s: graph validation failed: %s", path_prefix, e$message))
          })
        }
      }

      # Check pathway scoring
      if (is.null(pathway$scoring)) {
        warnings <- c(warnings, sprintf("%s: missing 'scoring' field", path_prefix))
      } else {
        if (!is.null(pathway$scoring$marker_mode)) {
          if (!pathway$scoring$marker_mode %in% c("any", "all")) {
            errors <- c(errors, sprintf("%s: marker_mode must be 'any' or 'all'", path_prefix))
          }
          # Check for marker genes
          marker_count <- sum(sapply(pathway$nodes, function(n) !is.null(n$marker) && n$marker == TRUE))
          if (marker_count == 0) {
            warnings <- c(warnings, sprintf("%s: marker_mode specified but no marker genes", path_prefix))
          }
        }
      }
    }
  }

  # Check that edges field is absent (edges belong in pathways)
  if (!is.null(data$edges) && length(data$edges) > 0) {
    errors <- c(errors, "Multi-pathway network should not have global 'edges' field (edges belong in individual pathways)")
  }

  # Strict checks
  if (strict) {
    if (is.null(data$source) || nchar(data$source) == 0) {
      warnings <- c(warnings, "No source specified")
    }
    if (is.null(data$tags) || length(data$tags) == 0) {
      warnings <- c(warnings, "No tags specified")
    }
  }

  list(valid = length(errors) == 0, errors = errors, warnings = warnings)
}

#' Print validation results
#'
#' @param validation_result Result from validate_potato()
#' @export
print_validation <- function(validation_result) {
  if (validation_result$valid) {
    cat("[VALID] Potato is valid\n")
  } else {
    cat("[ERROR] Potato has errors:\n")
    for (error in validation_result$errors) {
      cat("  ERROR:", error, "\n")
    }
  }

  if (length(validation_result$warnings) > 0) {
    cat("\nWarnings:\n")
    for (warning in validation_result$warnings) {
      cat("  WARNING:", warning, "\n")
    }
  }

  invisible(validation_result)
}


#' Print method for Potato
#' @export
S7::method(print, Potato) <- function(x, ...) {
  cat("<Potato:", x@id, ">\n")
  cat("  Name:", x@name, "\n")
  cat("  Nodes:", length(x@nodes), "\n")
  cat("  Edges:", length(x@edges), "\n")
  if (length(x@tags) > 0) {
    cat("  Tags:", paste(x@tags, collapse = ", "), "\n")
  }
  if (nchar(x@source) > 0) {
    cat("  Source:", x@source, "\n")
  }
  invisible(x)
}


#' Summary method for Potato
#' @export
S7::method(summary, Potato) <- function(object, ...) {
  cat("Potato:", object@name, "\n")
  cat("  ID:", object@id, "\n")
  cat("  Nodes:", length(object@nodes), "\n")
  cat("  Edges:", length(object@edges), "\n")

  if (length(object@tags) > 0) {
    cat("  Tags:", paste(object@tags, collapse = ", "), "\n")
  }

  if (nchar(object@source) > 0) {
    cat("  Source:", object@source, "\n")
  }

  # Show enzyme nodes
  enzyme_nodes <- get_enzyme_nodes(object)
  if (length(enzyme_nodes) > 0) {
    cat("\n  Enzyme nodes:\n")
    for (node in enzyme_nodes) {
      cat("    -", node$id, ":", node$name, "\n")
    }
  }

  # Validation
  validation <- validate_potato(object)
  if (validation$valid) {
    cat("\n  [VALID] Potato structure is valid\n")
  } else {
    cat("\n  [ERROR] Potato has validation errors:\n")
    for (error in validation$errors) {
      cat("      -", error, "\n")
    }
  }

  invisible(object)
}
