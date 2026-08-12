#' Print potato structure as text
#'
#' Shows pathway flow with genes, compounds, and structure
#'
#' @param potato Potato S7 object, v2 potato object, or path to JSON
#' @param compact Show compact one-line view (default: TRUE)
#' @param show_compounds Include compound names in flow (default: TRUE, ignored if compact)
#' @param show_databases Show database annotations (kofam, blast, hmm, etc.) (default: FALSE)
#' @param show_ec Show EC numbers (default: FALSE)
#'
#' @export
print_potato <- function(potato, compact = TRUE, show_compounds = TRUE, show_databases = FALSE, show_ec = FALSE) {

  # Load if path provided
  if (is.character(potato)) {
    # Check schema version to determine loader
    potato_data <- jsonlite::read_json(potato, simplifyVector = FALSE)
    if (!is.null(potato_data$schema_version) && potato_data$schema_version == "v2") {
      potato <- load_potato_v2(potato)
    } else {
      potato <- load_potato(potato)
    }
  }

  # Check if v2 schema
  is_v2 <- S7::S7_inherits(potato, PotatoV2)

  if (is_v2) {
    # Handle v2 schema inline
    cli::cli_h1(potato@name)
    cli::cli_text("Source: {.field {potato@source}}")
    cli::cli_text("Tags: {.val {paste(potato@tags, collapse = ', ')}}")
    if (!is.null(potato@notes) && nchar(potato@notes) > 0) {
      cli::cli_text("Notes: {.field {potato@notes}}")
    }
    cli::cli_rule()

    # Show genes
    cli::cli_h2("Genes ({length(potato@genes)})")
    for (gene in potato@genes) {
      cli::cli_text("{gene$name} [{.field {gene$id}}]")

      if (show_ec && !is.null(gene$ec) && length(gene$ec) > 0) {
        cli::cli_text("  {.dim EC: {paste(gene$ec, collapse = ', ')}}")
      }

      if (show_databases && !is.null(gene$databases)) {
        db_info <- character()
        if (!is.null(gene$databases$kofam)) {
          ko_terms <- unlist(gene$databases$kofam)
          db_info <- c(db_info, paste0("KOfam: ", paste(ko_terms, collapse = ", ")))
        }
        if (!is.null(gene$databases$blast)) {
          blast_terms <- unlist(gene$databases$blast)
          db_info <- c(db_info, paste0("BLAST: ", paste(blast_terms, collapse = ", ")))
        }
        if (!is.null(gene$databases$hmm)) {
          hmm_terms <- unlist(gene$databases$hmm)
          db_info <- c(db_info, paste0("HMM: ", paste(hmm_terms, collapse = ", ")))
        }

        if (length(db_info) > 0) {
          cli::cli_text("  {.dim ({paste(db_info, collapse = ', ')})}")
        }
      }

      if (!is.null(gene$reactions) && length(gene$reactions) > 0) {
        cli::cli_text("  {.dim Reactions: {paste(gene$reactions, collapse = ', ')}}")
      }
    }

    cli::cli_text("")

    # Show compounds
    if (!is.null(potato@compounds) && length(potato@compounds) > 0) {
      cli::cli_h2("Compounds ({length(potato@compounds)})")
      for (compound in potato@compounds) {
        cli::cli_text("{compound$name} [{.field {compound$id}}]")
      }
      cli::cli_text("")
    }

    # Show pathway(s)
    if (!is.null(potato@pathways)) {
      pathway_names <- names(potato@pathways)
      cli::cli_h2("Pathway{?s} ({length(pathway_names)})")

      for (pathway_id in pathway_names) {
        pathway <- potato@pathways[[pathway_id]]
        cli::cli_text("{.strong {pathway$name %||% pathway_id}}")

        if (!is.null(pathway$notes) && nchar(pathway$notes) > 0) {
          cli::cli_text("  {.dim {pathway$notes}}")
        }

        # Show edges
        if (!is.null(pathway$edges) && length(pathway$edges) > 0) {
          cli::cli_text("  {.dim {length(pathway$edges)} edges}")
          for (edge in head(pathway$edges, 3)) {
            from_label <- if (grepl("^C\\d+", edge$from)) {
              compound <- Find(function(c) c$id == edge$from, potato@compounds)
              if (!is.null(compound)) compound$name else edge$from
            } else {
              edge$from
            }
            to_label <- if (grepl("^C\\d+", edge$to)) {
              compound <- Find(function(c) c$id == edge$to, potato@compounds)
              if (!is.null(compound)) compound$name else edge$to
            } else {
              edge$to
            }
            rxn_str <- if (!is.null(edge$reaction)) paste0(" [", edge$reaction, "]") else ""
            cli::cli_text("    {from_label} → {to_label}{rxn_str}")
          }
          if (length(pathway$edges) > 3) {
            cli::cli_text("    {.dim ... and {length(pathway$edges) - 3} more}")
          }
        }

        cli::cli_text("")
      }
    }

    return(invisible(potato))
  }

  # Check if multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  if (is_network) {
    print_multi_pathway_network(potato, show_databases = show_databases, show_ec = show_ec, show_compounds = show_compounds)
    return(invisible(potato))
  }

  if (compact) {
    print_pathway_compact(potato, show_databases = show_databases, show_ec = show_ec)
    return(invisible(potato))
  }

  cli::cli_h1(potato@name)
  cli::cli_text("Source: {.field {potato@source}}")
  cli::cli_text("Tags: {.val {paste(potato@tags, collapse = ', ')}}")
  cli::cli_rule()

  # Group nodes by step
  steps <- list()
  for (node in potato@genes) {
    step_num <- node$step
    if (is.list(step_num)) step_num <- step_num[[1]]  # Handle bifunctional

    if (is.null(steps[[as.character(step_num)]])) {
      steps[[as.character(step_num)]] <- list()
    }
    steps[[as.character(step_num)]][[length(steps[[as.character(step_num)]]) + 1]] <- node
  }

  # Sort by step number
  step_numbers <- as.integer(names(steps))
  steps <- steps[order(step_numbers)]

  # Show input (if field exists and is populated)
  input_val <- tryCatch(potato@input, error = function(e) NULL)
  if (!is.null(input_val) && length(input_val) > 0) {
    cli::cli_alert_info("INPUT: {.emph {input_val$compound}} ({input_val$kegg_compound})")
    cli::cli_text("")
  }

  # Show each step
  for (step_name in names(steps)) {
    step_nodes <- steps[[step_name]]

    # Build step header
    markers <- sapply(step_nodes, function(n) n$marker %||% FALSE)
    required <- sapply(step_nodes, function(n) n$required %||% FALSE)

    step_label <- paste0("Step ", step_name)
    if (any(markers)) step_label <- paste0(step_label, " {.emph [MARKER]}")
    if (all(required)) step_label <- paste0(step_label, " {.emph [required]}")

    cli::cli_h2(step_label)

    # Show alternatives
    if (length(step_nodes) == 1) {
      # Single gene
      node <- step_nodes[[1]]
      format_gene(node, show_databases, show_ec)
    } else {
      # Multiple nodes = OR alternatives
      cli::cli_text("{.strong ALTERNATIVES} (any one):")
      for (node in step_nodes) {
        format_gene(node, show_databases, show_ec, indent = TRUE)
      }
    }

    # Show compound produced (from edges)
    if (show_compounds) {
      compounds <- get_step_output_compounds(potato, step_name)
      if (length(compounds) > 0) {
        for (cmp in compounds) {
          cli::cli_text("   {cli::symbol$arrow_down} {.field {cmp}}")
        }
      }
    }

    cli::cli_text("")
  }

  # Show output (if field exists and is populated)
  output_val <- tryCatch(potato@output, error = function(e) NULL)
  if (!is.null(output_val) && length(output_val) > 0) {
    cli::cli_alert_info("OUTPUT: {.emph {output_val$compound}} ({output_val$kegg_compound})")
  }

  # Show scoring
  cli::cli_rule("Scoring")
  cli::cli_text("Min fraction: {.val {potato@scoring$min_fraction %||% 0.75}}")
  cli::cli_text("Marker mode: {.val {potato@scoring$marker_mode %||% 'any'}}")

  if (!is.null(potato@notes) && nchar(potato@notes) > 0) {
    cli::cli_rule("Notes")
    cli::cli_text(potato@notes)
  }

  invisible(potato)
}


#' Print multi-pathway network (internal)
#' @noRd
print_multi_pathway_network <- function(potato, show_databases = FALSE, show_ec = FALSE, show_compounds = TRUE) {

  cli::cli_h1(potato@name)
  cli::cli_text("Source: {.field {potato@source}}")
  cli::cli_text("Tags: {.val {paste(potato@tags, collapse = ', ')}}")
  if (!is.null(potato@notes) && nchar(potato@notes) > 0) {
    cli::cli_text("Notes: {.field {potato@notes}}")
  }
  cli::cli_rule()

  # Extract pathway types from edges
  pathway_names <- names(potato@edges)

  cli::cli_alert_info("Multi-pathway network with {length(pathway_names)} pathway(s)")
  cli::cli_text("")

  # Print each pathway
  for (pathway_name in pathway_names) {
    pathway_edges <- potato@edges[[pathway_name]]

    # Check if edges exist and have type field
    if (length(pathway_edges) > 0 && !is.null(pathway_edges[[1]]$type)) {
      pathway_type <- pathway_edges[[1]]$type
      type_label <- ifelse(pathway_type == "variant", "[VARIANT]", "[INDEPENDENT]")

      cli::cli_h2("{pathway_name} {.emph {type_label}}")

      # Get unique nodes involved in this pathway
      node_ids <- unique(c(
        sapply(pathway_edges, function(e) e$from),
        sapply(pathway_edges, function(e) e$to)
      ))

      # Remove compound nodes (start with C)
      gene_ids <- node_ids[!grepl("^C\\d+", node_ids)]

      # Get pathway info
      pathway_info <- potato@pathway_info[[pathway_name]]

      # Show genes
      for (gene_id in gene_ids) {
        # Find node
        node <- NULL
        for (n in potato@genes) {
          if (n$id == gene_id) {
            node <- n
            break
          }
        }

        if (!is.null(node)) {
          # Get pathway-specific attributes
          is_marker <- !is.null(pathway_info$markers) && gene_id %in% pathway_info$markers
          is_required <- !is.null(pathway_info$required) && gene_id %in% pathway_info$required

          # Format gene display
          marker_str <- if (is_marker) "*" else ""
          required_str <- if (!is_required) "^" else ""

          cli::cli_text("  {node$name}{marker_str}{required_str} [{.field {node$id}}]")

          if (show_ec && !is.null(node$ec)) {
            cli::cli_text("    {.dim EC: {paste(node$ec, collapse = ', ')}}")
          }

          if (show_databases) {
            db_info <- character()
            if (!is.null(node$databases$kofam)) {
              ko_terms <- unlist(node$databases$kofam)
              db_info <- c(db_info, paste0("KOfam: ", paste(ko_terms, collapse = ", ")))
            }
            if (!is.null(node$databases$blast)) {
              blast_terms <- unlist(node$databases$blast)
              db_info <- c(db_info, paste0("BLAST: ", paste(blast_terms, collapse = ", ")))
            }
            if (!is.null(node$databases$hmm)) {
              hmm_terms <- unlist(node$databases$hmm)
              db_info <- c(db_info, paste0("HMM: ", paste(hmm_terms, collapse = ", ")))
            }

            if (length(db_info) > 0) {
              cli::cli_text("    {.dim ({paste(db_info, collapse = ', ')})}")
            }
          }
        }
      }

      # Show pathway flow if compounds requested
      if (show_compounds && length(pathway_edges) > 0) {
        cli::cli_text("")
        cli::cli_text("  {.strong Flow:}")
        for (edge in pathway_edges) {
          from_label <- if (grepl("^C\\d+", edge$from)) {
            paste0("<", edge$compound %||% edge$from, ">")
          } else {
            edge$from
          }
          to_label <- if (grepl("^C\\d+", edge$to)) {
            paste0("<", edge$compound %||% edge$to, ">")
          } else {
            edge$to
          }
          cli::cli_text("    {from_label} → {to_label}")
        }
      }

      cli::cli_text("")
    }
  }

  # Show legend
  cli::cli_rule("Legend")
  cli::cli_text("* = marker gene")
  cli::cli_text("^ = optional gene")

  invisible(potato)
}


#' Print pathway compact (internal)
#' @noRd
print_pathway_compact <- function(potato, show_databases = FALSE, show_ec = TRUE) {

  cli::cli_h1(potato@name)
  cli::cli_text("Source: {.field {potato@source}}")

  if (!is.null(potato@notes) && nchar(potato@notes) > 0) {
    cli::cli_text("Notes: {.field {potato@notes}}")
  }

  cli::cli_text("")

  # Show input
  input_val <- tryCatch(potato@input, error = function(e) NULL)
  if (!is.null(input_val) && length(input_val) > 0) {
    cli::cli_text("{.emph <{input_val$compound} [{input_val$kegg_compound}]>} → ", appendLF = FALSE)
  }

  # Group nodes by step
  steps <- list()
  for (node in potato@genes) {
    step_num <- node$step
    if (is.list(step_num)) step_num <- step_num[[1]]

    if (is.null(steps[[as.character(step_num)]])) {
      steps[[as.character(step_num)]] <- list()
    }
    steps[[as.character(step_num)]][[length(steps[[as.character(step_num)]]) + 1]] <- node
  }

  # Sort by step
  step_numbers <- as.integer(names(steps))
  steps <- steps[order(step_numbers)]

  # Show flow
  flow_parts <- character()
  for (step_name in names(steps)) {
    step_nodes <- steps[[step_name]]

    if (length(step_nodes) == 1) {
      flow_parts <- c(flow_parts, format_gene_compact(step_nodes[[1]], show_databases, show_ec))
    } else {
      # OR alternatives
      alts <- sapply(step_nodes, function(n) format_gene_compact(n, show_databases, show_ec))
      flow_parts <- c(flow_parts, paste0("(", paste(alts, collapse = " | "), ")"))
    }

    # Add compound if available
    compounds <- get_step_output_compounds(potato, step_name)
    if (length(compounds) > 0) {
      for (cmp in compounds) {
        flow_parts <- c(flow_parts, paste0("<", cmp, ">"))
      }
    }
  }

  cli::cli_text(paste(flow_parts, collapse = " → "))

  # Show output
  output_val <- tryCatch(potato@output, error = function(e) NULL)
  if (!is.null(output_val) && length(output_val) > 0) {
    cli::cli_text(" → {.emph <{output_val$compound} [{output_val$kegg_compound}]>}")
  } else {
    cli::cli_text("")
  }

  cli::cli_text("")
  cli::cli_text("{.dim * = marker, ^ = optional}")

  invisible(potato)
}


#' Format gene compact (internal)
#' @noRd
format_gene_compact <- function(node, show_databases = FALSE, show_ec = TRUE) {
  marker_str <- if (node$marker %||% FALSE) "*" else ""
  optional_str <- if (!(node$required %||% TRUE)) "^" else ""

  label <- paste0(node$id, marker_str, optional_str)

  if (show_ec && !is.null(node$ec) && length(node$ec) > 0) {
    label <- paste0(label, "[", paste(node$ec, collapse = ","), "]")
  }

  if (show_databases) {
    db_parts <- character()
    if (!is.null(node$databases$kofam)) {
      db_parts <- c(db_parts, paste0("K:", paste(unlist(node$databases$kofam), collapse = ",")))
    }
    if (!is.null(node$databases$blast)) {
      db_parts <- c(db_parts, paste0("B:", length(node$databases$blast)))
    }
    if (!is.null(node$databases$hmm)) {
      db_parts <- c(db_parts, paste0("H:", length(node$databases$hmm)))
    }

    if (length(db_parts) > 0) {
      label <- paste0(label, "{", paste(db_parts, collapse = ","), "}")
    }
  }

  label
}


#' Format gene for printing (internal)
#' @noRd
format_gene <- function(node, show_databases, show_ec, indent = FALSE) {
  marker_str <- if (node$marker %||% FALSE) " {.emph [MARKER]}" else ""
  optional_str <- if (!(node$required %||% TRUE)) " {.emph [optional]}" else ""

  spaces <- if (indent) "  " else ""

  cli::cli_text(paste0(spaces, "{node$name} [{.field {node$id}}]", marker_str, optional_str))

  # Show EC
  if (show_ec && !is.null(node$ec)) {
    spaces2 <- if (indent) "    " else "  "
    cli::cli_text(paste0(spaces2, "{.dim EC: {paste(node$ec, collapse = ', ')}}"))
  }

  # Show databases
  if (show_databases) {
    db_info <- character()
    if (!is.null(node$databases$kofam)) {
      ko_terms <- unlist(node$databases$kofam)
      db_info <- c(db_info, paste0("KOfam: ", paste(ko_terms, collapse = ", ")))
    }
    if (!is.null(node$databases$blast)) {
      blast_terms <- unlist(node$databases$blast)
      db_info <- c(db_info, paste0("BLAST: ", paste(blast_terms, collapse = ", ")))
    }
    if (!is.null(node$databases$hmm)) {
      hmm_terms <- unlist(node$databases$hmm)
      db_info <- c(db_info, paste0("HMM: ", paste(hmm_terms, collapse = ", ")))
    }

    if (length(db_info) > 0) {
      spaces <- if (indent) "      " else "    "
      cli::cli_text(paste0(spaces, "{.dim (", paste(db_info, collapse = ", "), ")}"))
    }
  }
}


#' Get output compounds for a step (internal)
#' @noRd
get_step_output_compounds <- function(potato, step_name) {
  step_num <- as.integer(step_name)
  compounds <- character()

  for (edge in potato@edges) {
    from_step <- as.integer(sub(".*_(\\d+)$", "\\1", edge$from))
    if (from_step == step_num && !is.null(edge$compound)) {
      if (!(edge$compound %in% compounds)) {
        compounds <- c(compounds, edge$compound)
      }
    }
  }

  compounds
}


