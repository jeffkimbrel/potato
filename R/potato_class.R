#' Potato S7 Class
#'
#' Represents a pathway definition with genes and pathway structure
#'
#' @param id Character. Unique identifier for the potato
#' @param name Character. Human-readable name
#' @param type Character. Type of biological entity. Options:
#'   - "pathway" (default): Metabolic pathway
#'   - "structural_complex": Physical structures (R-bodies, flagella, pili)
#'   - "secretion_system": Type III/IV/VI secretion systems
#'   - "regulatory_system": Two-component systems, quorum sensing
#'   - "defense_system": CRISPR, restriction-modification
#'   - "transport_system": ABC transporters, porins
#'   - "marker_gene": Single diagnostic genes
#' @param nodes List. Nodes (genes/enzymes) in the pathway. Each node can have:
#'   - `marker`: Logical. If TRUE, this gene is a diagnostic marker for the pathway
#'   - `required`: Logical. If TRUE, pathway cannot be complete without this gene
#'   - `type`: Character. "enzyme" or "compound"
#'   - Detection methods: `ko`, `ec`, `pfam`, `hmm`, `blast_db`, `blast_terms`
#'   - Direct path overrides: `blast_db_path`, `hmm_path`, `pfam_path` (override config lookup)
#'   - `thresholds`: List. Gene-specific thresholds overriding config defaults
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
    type = S7::class_character,
    nodes = S7::class_list,
    edges = S7::class_list,
    tags = S7::class_character,
    source = S7::class_character,
    notes = S7::class_character,
    scoring = S7::class_list,
    json_path = S7::class_character,
    graph = S7::class_any
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
#' @param path Path to potato JSON file
#' @return Potato object
#' @export
load_potato <- function(path) {
  if (!file.exists(path)) {
    stop("Potato file not found: ", path, call. = FALSE)
  }

  # Don't simplify vectors - keep list structure
  data <- jsonlite::read_json(path, simplifyVector = FALSE)

  # Convert nodes to list of lists (keep structure)
  nodes <- if (is.null(data$nodes)) list() else data$nodes

  # Convert edges to list of lists
  edges <- if (is.null(data$edges)) list() else data$edges

  # Flatten simple vectors
  tags <- if (is.null(data$tags)) character(0) else unlist(data$tags)

  Potato(
    id = data$id,
    name = data$name,
    type = if (is.null(data$type)) "pathway" else data$type,
    nodes = nodes,
    edges = edges,
    tags = tags,
    source = if (is.null(data$source)) "" else data$source,
    notes = if (is.null(data$notes)) "" else data$notes,
    scoring = if (is.null(data$scoring)) list() else data$scoring,
    json_path = path,
    graph = NULL  # Build on demand
  )
}


#' Load multiple potatoes from directory
#'
#' @param dir Directory containing potato JSON files
#' @param tags Optional character vector of tags to filter by
#' @return Named list of Potato objects
#' @export
load_potatoes <- function(dir, tags = NULL) {
  if (!dir.exists(dir)) {
    stop("Directory not found: ", dir, call. = FALSE)
  }

  json_files <- list.files(dir, pattern = "\\.json$", full.names = TRUE)

  potatoes <- lapply(json_files, function(filepath) {
    tryCatch({
      load_potato(filepath)
    }, error = function(e) {
      warning("Failed to load ", basename(filepath), ": ", e$message)
      NULL
    })
  })

  # Remove NULLs from failed loads
  potatoes <- potatoes[!sapply(potatoes, is.null)]

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
#' @export
get_enzyme_nodes <- function(potato) {
  S7::check_is_S7(potato)

  # Filter nodes that are enzymes
  enzyme_nodes <- Filter(function(node) {
    !is.null(node$type) && node$type == "enzyme"
  }, potato@nodes)

  enzyme_nodes
}


#' Get detection terms for a database
#'
#' Extracts all terms (KO IDs, gene names, etc.) for a specific database from
#' a potato's enzyme nodes. Supports both new 'databases' field and legacy fields.
#'
#' @param potato Potato object
#' @param database_name Name of database in config (e.g., "kofam113", "gator_blast")
#' @param tool_type Legacy parameter for backward compatibility ("ko", "pfam", "ec", "hmm", "blast_db", "blast_terms")
#' @return Character vector of unique terms
#' @export
get_detection_terms <- function(potato, database_name = NULL, tool_type = NULL) {
  S7::check_is_S7(potato)

  if (is.null(database_name) && is.null(tool_type)) {
    stop("Must provide either database_name or tool_type", call. = FALSE)
  }

  enzyme_nodes <- get_enzyme_nodes(potato)

  # New schema: extract terms from databases field
  if (!is.null(database_name)) {
    terms <- unlist(lapply(enzyme_nodes, function(node) {
      if (!is.null(node$databases) && !is.null(node$databases[[database_name]])) {
        node$databases[[database_name]]
      }
    }))
    return(unique(terms))
  }

  # Legacy schema: extract from tool_type field
  if (!is.null(tool_type)) {
    terms <- unlist(lapply(enzyme_nodes, function(node) {
      if (!is.null(node[[tool_type]])) {
        node[[tool_type]]
      }
    }))
    return(unique(terms))
  }
}


#' Get marker gene nodes from potato
#'
#' @param potato Potato object
#' @return List of nodes marked as marker genes
#' @export
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
#' @export
build_potato_graph <- function(potato) {
  S7::check_is_S7(potato)

  if (length(potato@edges) == 0) {
    # Graph with just nodes, no edges
    node_ids <- sapply(potato@nodes, function(n) n$id)
    g <- igraph::make_empty_graph(n = length(node_ids), directed = TRUE)
    igraph::V(g)$name <- node_ids
    return(g)
  }

  # Build edge list
  edge_list <- do.call(rbind, lapply(potato@edges, function(e) {
    c(e$from, e$to)
  }))

  # Create graph
  g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)

  # That's it - just return the basic graph structure
  # Can add attributes later if needed

  g
}


#' Validate potato structure
#'
#' Validates a potato JSON structure for common errors and required fields.
#' Users should run this on custom potatoes before using them for annotation.
#'
#' @param potato Potato object
#' @param strict Logical. If TRUE, performs additional strict checks (default FALSE)
#' @return List with 'valid' (logical), 'errors' (character vector), and 'warnings' (character vector)
#' @export
validate_potato <- function(potato, strict = FALSE) {
  S7::check_is_S7(potato)

  errors <- character(0)
  warnings <- character(0)

  # Check required fields
  if (nchar(potato@id) == 0) {
    errors <- c(errors, "Missing required field: 'id'")
  }
  if (nchar(potato@name) == 0) {
    errors <- c(errors, "Missing required field: 'name'")
  }

  # Validate type field
  valid_types <- c("pathway", "structural_complex", "secretion_system",
                   "regulatory_system", "defense_system", "transport_system",
                   "marker_gene")
  if (nchar(potato@type) > 0 && !potato@type %in% valid_types) {
    warnings <- c(warnings, sprintf("Unknown type '%s' (expected: %s)",
                                   potato@type, paste(valid_types, collapse = ", ")))
  }

  # Check ID format (alphanumeric, underscores, hyphens only)
  if (!grepl("^[a-zA-Z0-9_-]+$", potato@id)) {
    warnings <- c(warnings, "Potato ID should contain only letters, numbers, underscores, and hyphens")
  }

  # Check nodes
  if (length(potato@nodes) == 0) {
    errors <- c(errors, "Potato has no nodes")
  }

  # Check node IDs are unique
  node_ids <- sapply(potato@nodes, function(n) n$id)
  if (length(node_ids) != length(unique(node_ids))) {
    errors <- c(errors, "Duplicate node IDs found")
  }

  # Validate each node
  for (i in seq_along(potato@nodes)) {
    node <- potato@nodes[[i]]
    node_prefix <- sprintf("Node %d (%s)", i, node$id %||% "unnamed")

    # Required node fields
    if (is.null(node$id)) {
      errors <- c(errors, sprintf("%s: missing 'id'", node_prefix))
    }
    if (is.null(node$type)) {
      errors <- c(errors, sprintf("%s: missing 'type'", node_prefix))
    } else if (!node$type %in% c("enzyme", "compound")) {
      warnings <- c(warnings, sprintf("%s: type '%s' is non-standard (expected 'enzyme' or 'compound')",
                                      node_prefix, node$type))
    }

    # Check enzyme nodes have detection methods via databases field
    if (!is.null(node$type) && node$type == "enzyme") {
      # New schema: check for databases field
      has_databases <- !is.null(node$databases) && length(node$databases) > 0
      # Legacy schema: check for old detection method fields
      has_legacy <- any(c("ko", "ec", "pfam", "hmm", "blast_db", "blast_terms") %in% names(node))

      if (!has_databases && !has_legacy) {
        warnings <- c(warnings, sprintf("%s: enzyme has no detection methods (needs 'databases' field or legacy fields)",
                                        node_prefix))
      }

      # Validate databases field structure
      if (has_databases) {
        if (!is.list(node$databases)) {
          errors <- c(errors, sprintf("%s: 'databases' must be a list", node_prefix))
        } else {
          # Check each database entry
          for (db_name in names(node$databases)) {
            db_terms <- node$databases[[db_name]]
            if (!is.character(db_terms) && !is.list(db_terms)) {
              errors <- c(errors, sprintf("%s: database '%s' terms must be character vector or list",
                                          node_prefix, db_name))
            }
            # Validate KO format if this looks like a kofam database
            if (grepl("kofam", db_name, ignore.case = TRUE)) {
              ko_ids <- if (is.list(db_terms)) unlist(db_terms) else db_terms
              invalid_ko <- ko_ids[!grepl("^K[0-9]{5}$", ko_ids)]
              if (length(invalid_ko) > 0) {
                warnings <- c(warnings, sprintf("%s: invalid KO format in '%s': %s (should be K##### e.g., K00001)",
                                                node_prefix, db_name, paste(invalid_ko, collapse = ", ")))
              }
            }
          }
        }
      }

      # Legacy field validation (keep for backward compatibility)
      if (!is.null(node$ko)) {
        ko_ids <- if (is.list(node$ko)) unlist(node$ko) else node$ko
        invalid_ko <- ko_ids[!grepl("^K[0-9]{5}$", ko_ids)]
        if (length(invalid_ko) > 0) {
          warnings <- c(warnings, sprintf("%s: invalid KO format: %s (should be K##### e.g., K00001)",
                                          node_prefix, paste(invalid_ko, collapse = ", ")))
        }
      }

      if (!is.null(node$ec)) {
        ec_ids <- if (is.list(node$ec)) unlist(node$ec) else node$ec
        invalid_ec <- ec_ids[!grepl("^[0-9]+\\.[0-9-]+\\.[0-9-]+\\.[0-9-]+$", ec_ids)]
        if (length(invalid_ec) > 0) {
          warnings <- c(warnings, sprintf("%s: invalid EC format: %s (should be X.X.X.X)",
                                          node_prefix, paste(invalid_ec, collapse = ", ")))
        }
      }

      # Check required field type
      if (!is.null(node$required) && !is.logical(node$required)) {
        errors <- c(errors, sprintf("%s: 'required' must be TRUE/FALSE, not '%s'",
                                    node_prefix, node$required))
      }

      # Check marker field type
      if (!is.null(node$marker) && !is.logical(node$marker)) {
        errors <- c(errors, sprintf("%s: 'marker' must be TRUE/FALSE, not '%s'",
                                    node_prefix, node$marker))
      }

      # Validate thresholds
      if (!is.null(node$thresholds)) {
        valid_threshold_names <- c("kofam_score", "kofam_evalue", "blast_evalue",
                                   "blast_bitscore", "blast_pident", "hmm_evalue",
                                   "hmm_bitscore", "hmm_domain_evalue")
        invalid_thresholds <- setdiff(names(node$thresholds), valid_threshold_names)
        if (length(invalid_thresholds) > 0) {
          warnings <- c(warnings, sprintf("%s: unknown threshold names: %s",
                                          node_prefix, paste(invalid_thresholds, collapse = ", ")))
        }
      }

      # Check for direct path overrides (legacy support)
      if (!is.null(node$blast_db_path)) {
        if (!file.exists(node$blast_db_path)) {
          warnings <- c(warnings, sprintf("%s: blast_db_path specified but file not found: %s",
                                          node_prefix, node$blast_db_path))
        }
        if (has_databases) {
          warnings <- c(warnings, sprintf("%s: both 'databases' and 'blast_db_path' specified (path takes precedence)",
                                          node_prefix))
        }
      }

      if (!is.null(node$hmm_path)) {
        if (!file.exists(node$hmm_path)) {
          warnings <- c(warnings, sprintf("%s: hmm_path specified but file not found: %s",
                                          node_prefix, node$hmm_path))
        }
        if (has_databases) {
          warnings <- c(warnings, sprintf("%s: both 'databases' and 'hmm_path' specified (path takes precedence)",
                                          node_prefix))
        }
      }

      if (!is.null(node$pfam_path)) {
        if (!file.exists(node$pfam_path)) {
          warnings <- c(warnings, sprintf("%s: pfam_path specified but file not found: %s",
                                          node_prefix, node$pfam_path))
        }
        if (has_databases) {
          warnings <- c(warnings, sprintf("%s: both 'databases' and 'pfam_path' specified (path takes precedence)",
                                          node_prefix))
        }
      }
    }
  }

  # Check edges reference valid nodes
  for (i in seq_along(potato@edges)) {
    edge <- potato@edges[[i]]

    if (is.null(edge$from)) {
      errors <- c(errors, sprintf("Edge %d: missing 'from' field", i))
    } else if (!edge$from %in% node_ids) {
      errors <- c(errors, sprintf("Edge %d: 'from' references non-existent node '%s'", i, edge$from))
    }

    if (is.null(edge$to)) {
      errors <- c(errors, sprintf("Edge %d: missing 'to' field", i))
    } else if (!edge$to %in% node_ids) {
      errors <- c(errors, sprintf("Edge %d: 'to' references non-existent node '%s'", i, edge$to))
    }
  }

  # Check for cycles (DAG requirement)
  if (length(errors) == 0 && length(potato@edges) > 0) {
    tryCatch({
      g <- build_potato_graph(potato)
      if (!igraph::is_dag(g)) {
        errors <- c(errors, "Potato contains cycles (must be a directed acyclic graph)")
      }
    }, error = function(e) {
      errors <- c(errors, sprintf("Graph construction failed: %s", e$message))
    })
  }

  # Check for at least one marker gene if marker_mode is specified
  if (!is.null(potato@scoring$marker_mode)) {
    marker_nodes <- Filter(function(n) !is.null(n$marker) && n$marker == TRUE, potato@nodes)
    if (length(marker_nodes) == 0) {
      warnings <- c(warnings, "Scoring has marker_mode specified but no nodes have marker=true")
    }

    # Validate marker_mode value
    if (!potato@scoring$marker_mode %in% c("any", "all")) {
      errors <- c(errors, "scoring.marker_mode must be 'any' or 'all'")
    }
  }

  # Strict mode checks
  if (strict) {
    # Check for source attribution
    if (nchar(potato@source) == 0) {
      warnings <- c(warnings, "No 'source' specified (good practice to cite pathway origin)")
    }

    # Check for tags
    if (length(potato@tags) == 0) {
      warnings <- c(warnings, "No tags specified (tags help organize potato collections)")
    }

    # Check for scoring parameters
    if (length(potato@scoring) == 0) {
      warnings <- c(warnings, "No scoring parameters specified")
    }

    # Recommend marker genes for complex pathways
    enzyme_count <- length(Filter(function(n) !is.null(n$type) && n$type == "enzyme", potato@nodes))
    marker_count <- length(Filter(function(n) !is.null(n$marker) && n$marker == TRUE, potato@nodes))
    if (enzyme_count >= 5 && marker_count == 0) {
      warnings <- c(warnings, "Complex pathway (5+ enzymes) without marker genes - consider designating key diagnostic genes")
    }
  }

  list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings
  )
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
