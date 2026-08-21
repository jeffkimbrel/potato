#' Potato S7 Class
#'
#' Represents a pathway definition with genes and pathway structure
#'
#' @param id Character. Unique identifier for the potato
#' @param name Character. Human-readable name
#' @param genes List. Genes (enzymes) in the pathway. Each gene can have:
#'   - `marker`: Logical. If TRUE, this gene is a diagnostic marker for the pathway
#'   - `required`: Logical. If TRUE, pathway cannot be complete without this gene
#'   - `type`: Character. "enzyme", "compound", or "transporter"
#'   - `databases`: List. Detection methods by database type (e.g., kofam, blast, hmm)
#' @param edges List. Edges connecting genes or pathways object for multi-pathway networks
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
    genes = S7::class_list,
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
    } else if (length(self@genes) == 0) {
      "Potato must have at least one gene"
    }
  }
)



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
    cli::cli_abort("Directory not found: {.path {dir}}")
  }

  json_files <- list.files(dir, pattern = "\\.json$", full.names = TRUE)

  potatoes <- lapply(json_files, function(filepath) {
    tryCatch({
      # Read JSON to check schema version and active flag
      data <- jsonlite::read_json(filepath, simplifyVector = FALSE)

      # Skip inactive potatoes unless requested
      if (!include_inactive && !is.null(data$active) && data$active == FALSE) {
        return(NULL)
      }

      return(load_potato_v2(filepath))
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

  # Name by ID (handle both v1 Potato and v2 PotatoV2 S7 objects)
  names(potatoes) <- sapply(potatoes, function(p) p@id)

  # Filter by tags if specified
  if (!is.null(tags)) {
    potatoes <- potatoes[sapply(potatoes, function(p) {
      any(tags %in% p@tags)
    })]
  }

  potatoes
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
  }, potato@genes)

  enzyme_nodes
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
      # Check if PotatoV2 (has schema_version) or old Potato class
      is_v2 <- !is.null(tryCatch(potato@schema_version, error = function(e) NULL))

      if (is_v2) {
        # PotatoV2 object
        data <- list(
          schema_version = potato@schema_version,
          id = potato@id,
          name = potato@name,
          genes = potato@genes,
          compounds = potato@compounds,
          pathways = potato@pathways,
          tags = potato@tags,
          source = potato@source,
          notes = potato@notes
        )
      } else {
        # Old Potato class (should not exist anymore, but handle gracefully)
        data <- list(
          id = potato@id,
          name = potato@name,
          genes = potato@genes,
          tags = potato@tags,
          source = potato@source
        )

        # Add edges/pathways if present
        if (!is.null(tryCatch(potato@edges, error = function(e) NULL))) {
          data$edges <- potato@edges
        }
      }
    } else {
      data <- potato
    }
  } else {
    cli::cli_abort("{.arg potato} must be a {.cls Potato} object or list")
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

  result <- validate_multi_pathway(data, strict)
  return(result)

}


#' Validate multi-pathway network potato (V2 schema)
#' @keywords internal
validate_multi_pathway <- function(data, strict) {
  errors <- character(0)
  warnings <- character(0)

  # Check genes array exists (not "nodes" - that was old schema)
  if (is.null(data$genes) || length(data$genes) == 0) {
    errors <- c(errors, "Multi-pathway network has no global genes")
  }

  # Collect gene IDs and compound IDs for edge validation
  gene_ids <- character()
  compound_ids <- character()

  if (!is.null(data$genes)) {
    gene_ids <- sapply(data$genes, function(g) g$id)

    # Check gene IDs are unique
    if (length(gene_ids) != length(unique(gene_ids))) {
      errors <- c(errors, "Duplicate gene IDs in global genes array")
    }

    # Validate global genes (detection methods only, no step/required/marker)
    for (i in seq_along(data$genes)) {
      gene <- data$genes[[i]]
      gene_prefix <- sprintf("Global gene %d (%s)", i, gene$id %||% "unnamed")

      if (is.null(gene$id)) {
        errors <- c(errors, sprintf("%s: missing 'id'", gene_prefix))
      }
      if (is.null(gene$type)) {
        warnings <- c(warnings, sprintf("%s: missing 'type' (recommended)", gene_prefix))
      }

      # Global genes should NOT have pathway-specific attributes
      if (!is.null(gene$step)) {
        errors <- c(errors, sprintf("%s: 'step' should not be in V2 genes (step-based scoring removed)", gene_prefix))
      }
      if (!is.null(gene$required)) {
        errors <- c(errors, sprintf("%s: 'required' belongs on edges, not genes", gene_prefix))
      }
      if (!is.null(gene$marker)) {
        errors <- c(errors, sprintf("%s: 'marker' belongs on edges, not genes", gene_prefix))
      }

      # Validate databases field
      if (!is.null(gene$type) && gene$type == "enzyme") {
        if (is.null(gene$databases) || length(gene$databases) == 0) {
          warnings <- c(warnings, sprintf("%s: enzyme missing 'databases' field (no detection methods)", gene_prefix))
        } else {
          valid_db_types <- c("kofam", "blast", "hmm")
          for (db_name in names(gene$databases)) {
            if (!db_name %in% valid_db_types) {
              errors <- c(errors, sprintf("%s: invalid database type '%s' (must be kofam, blast, or hmm)", gene_prefix, db_name))
            }

            # Validate KO format
            if (db_name == "kofam") {
              ko_ids <- unlist(gene$databases[[db_name]])
              invalid_ko <- ko_ids[!grepl("^K[0-9]{5}$", ko_ids)]
              if (length(invalid_ko) > 0) {
                warnings <- c(warnings, sprintf("%s: invalid KO format: %s", gene_prefix, paste(invalid_ko, collapse = ", ")))
              }
            }
          }
        }
      }
    }
  }

  # Collect compound IDs
  if (!is.null(data$compounds)) {
    compound_ids <- sapply(data$compounds, function(c) c$id)

    # Check compound IDs are unique
    if (length(compound_ids) != length(unique(compound_ids))) {
      errors <- c(errors, "Duplicate compound IDs in compounds array")
    }
  }

  # Valid node IDs = genes + compounds
  valid_node_ids <- c(gene_ids, compound_ids)

  # Validate pathways field
  if (is.null(data$pathways) || !is.list(data$pathways) || length(data$pathways) == 0) {
    errors <- c(errors, "Multi-pathway network must have 'pathways' field with at least one pathway")
  } else {
    # Validate each pathway
    for (pathway_id in names(data$pathways)) {
      pathway <- data$pathways[[pathway_id]]
      path_prefix <- sprintf("Pathway '%s'", pathway_id)

      # Check pathway has edges (not nodes - that was old schema)
      if (is.null(pathway$edges) || length(pathway$edges) == 0) {
        warnings <- c(warnings, sprintf("%s: no edges defined (empty pathway)", path_prefix))
      }

      # Check pathway type
      if (is.null(pathway$type)) {
        warnings <- c(warnings, sprintf("%s: missing 'type' field (should be 'variant' or 'independent')", path_prefix))
      } else if (!pathway$type %in% c("variant", "independent")) {
        errors <- c(errors, sprintf("%s: type must be 'variant' or 'independent', got '%s'", path_prefix, pathway$type))
      }

      # Validate pathway edges
      if (!is.null(pathway$edges) && length(pathway$edges) > 0) {
        for (i in seq_along(pathway$edges)) {
          edge <- pathway$edges[[i]]
          edge_prefix <- sprintf("%s edge %d", path_prefix, i)

          # Validate 'from' field
          from_is_null <- is.null(edge$from) || (is.list(edge$from) && length(edge$from) == 0)
          if (from_is_null) {
            errors <- c(errors, sprintf("%s: missing 'from'", edge_prefix))
          } else {
            if (!is.character(edge$from) || length(edge$from) != 1) {
              errors <- c(errors, sprintf("%s: 'from' must be a single string, got type %s", edge_prefix, class(edge$from)[1]))
            } else if (nchar(edge$from) == 0) {
              errors <- c(errors, sprintf("%s: 'from' is empty string", edge_prefix))
            } else if (!edge$from %in% valid_node_ids) {
              errors <- c(errors, sprintf("%s: 'from' node '%s' not found in genes or compounds", edge_prefix, edge$from))
            }
          }

          # Validate 'to' field
          to_is_null <- is.null(edge$to) || (is.list(edge$to) && length(edge$to) == 0)
          if (to_is_null) {
            errors <- c(errors, sprintf("%s: missing 'to'", edge_prefix))
          } else {
            if (!is.character(edge$to) || length(edge$to) != 1) {
              errors <- c(errors, sprintf("%s: 'to' must be a single string, got type %s", edge_prefix, class(edge$to)[1]))
            } else if (nchar(edge$to) == 0) {
              errors <- c(errors, sprintf("%s: 'to' is empty string", edge_prefix))
            } else if (!edge$to %in% valid_node_ids) {
              errors <- c(errors, sprintf("%s: 'to' node '%s' not found in genes or compounds", edge_prefix, edge$to))
            }
          }

          # Validate required/marker attributes on edges
          if (!is.null(edge$required) && !is.logical(edge$required)) {
            errors <- c(errors, sprintf("%s: 'required' must be TRUE/FALSE", edge_prefix))
          }
          if (!is.null(edge$marker) && !is.logical(edge$marker)) {
            errors <- c(errors, sprintf("%s: 'marker' must be TRUE/FALSE", edge_prefix))
          }

          # Validate reaction if present
          if (!is.null(edge$reaction)) {
            if (!is.character(edge$reaction) || length(edge$reaction) != 1) {
              errors <- c(errors, sprintf("%s: 'reaction' must be a single string (e.g., 'R00001')", edge_prefix))
            }
          }
        }

        # Check for cycles in this pathway
        if (length(errors) == 0) {
          tryCatch({
            edge_list <- do.call(rbind, lapply(pathway$edges, function(e) {
              if (!is.null(e$from) && !is.null(e$to)) c(e$from, e$to) else NULL
            }))

            if (!is.null(edge_list) && nrow(edge_list) > 0) {
              g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)
              if (!igraph::is_dag(g)) {
                warnings <- c(warnings, sprintf("%s: contains cycles (metabolic cycles are OK, but check carefully)", path_prefix))
              }
            }
          }, error = function(e) {
            warnings <- c(warnings, sprintf("%s: could not check for cycles: %s", path_prefix, e$message))
          })
        }
      }

      # Validate input/output kegg_compound
      if (!is.null(pathway$input) && !is.null(pathway$input$kegg_compound)) {
        kegg_id <- pathway$input$kegg_compound
        if (!is.character(kegg_id) || length(kegg_id) != 1) {
          errors <- c(errors, sprintf("%s input: 'kegg_compound' must be a single string (e.g., 'C00024')", path_prefix))
        } else if (!grepl("^C[0-9]{5}$", kegg_id)) {
          warnings <- c(warnings, sprintf("%s input: 'kegg_compound' format should be C##### (e.g., 'C00024'), got '%s'", path_prefix, kegg_id))
        }
      }

      if (!is.null(pathway$output) && !is.null(pathway$output$kegg_compound)) {
        kegg_id <- pathway$output$kegg_compound
        if (!is.character(kegg_id) || length(kegg_id) != 1) {
          errors <- c(errors, sprintf("%s output: 'kegg_compound' must be a single string (e.g., 'C00024')", path_prefix))
        } else if (!grepl("^C[0-9]{5}$", kegg_id)) {
          warnings <- c(warnings, sprintf("%s output: 'kegg_compound' format should be C##### (e.g., 'C00024'), got '%s'", path_prefix, kegg_id))
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
          # Check for marker edges
          if (!is.null(pathway$edges)) {
            marker_count <- sum(sapply(pathway$edges, function(e) !is.null(e$marker) && e$marker == TRUE))
            if (marker_count == 0) {
              warnings <- c(warnings, sprintf("%s: marker_mode specified but no marker edges", path_prefix))
            }
          }
        }

        # Check min_fraction against gene count for float precision traps
        if (!is.null(pathway$scoring$min_fraction) && !is.null(pathway$edges)) {
          min_frac <- pathway$scoring$min_fraction
          # Count unique genes in this pathway (exclude compounds)
          gene_nodes <- unique(unlist(lapply(pathway$edges, function(e) {
            c(if (!is.null(e$from) && e$from %in% gene_ids) e$from,
              if (!is.null(e$to) && e$to %in% gene_ids) e$to)
          })))
          n_genes <- length(gene_nodes)

          if (n_genes > 0) {
            # Calculate valid fractions for this gene count
            valid_fractions <- sapply(1:n_genes, function(i) i / n_genes)

            # Check if min_fraction falls between valid fractions
            for (i in seq_along(valid_fractions)[-length(valid_fractions)]) {
              lower <- valid_fractions[i]
              upper <- valid_fractions[i + 1]
              if (min_frac > lower && min_frac < upper) {
                warnings <- c(warnings, sprintf(
                  "%s: min_fraction %.2f requires %d/%d genes (%.1f%%) to pass. Did you mean %.2f for %d/%d genes (%.1f%%)?",
                  path_prefix, min_frac, i + 1, n_genes, upper * 100,
                  round(lower, 2), i, n_genes, lower * 100
                ))
              }
            }
          }
        }
      }
    }
  }

  # Check that edges field is absent at top level (edges belong in pathways)
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


#' Print method for Potato
#' @export
S7::method(print, Potato) <- function(x, ...) {
  cat("<Potato:", x@id, ">\n")
  cat("  Name:", x@name, "\n")
  cat("  Genes:", length(x@genes), "\n")
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
  cat("  Genes:", length(object@genes), "\n")
  cat("  Edges:", length(object@edges), "\n")

  if (length(object@tags) > 0) {
    cat("  Tags:", paste(object@tags, collapse = ", "), "\n")
  }

  if (nchar(object@source) > 0) {
    cat("  Source:", object@source, "\n")
  }

  # Show enzyme genes
  enzyme_genes <- get_enzyme_nodes(object)
  if (length(enzyme_genes) > 0) {
    cat("\n  Enzyme genes:\n")
    for (gene in enzyme_genes) {
      cat("    -", gene$id, ":", gene$name, "\n")
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
