#' Print potato structure as text
#'
#' Shows pathway flow with genes, compounds, and structure
#'
#' @param potato Potato S7 object or path to JSON
#' @param compact Show compact one-line view (default: TRUE)
#' @param show_compounds Include compound names in flow (default: TRUE, ignored if compact)
#' @param show_ko Show KO IDs (default: FALSE)
#' @param show_ec Show EC numbers (default: FALSE)
#'
#' @export
print_potato <- function(potato, compact = TRUE, show_compounds = TRUE, show_ko = FALSE, show_ec = FALSE) {

  # Load if path provided
  if (is.character(potato)) {
    potato <- load_potato(potato)
  }

  # Check if multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  if (is_network) {
    print_multi_pathway_network(potato, show_ko = show_ko, show_ec = show_ec)
    return(invisible(potato))
  }

  if (compact) {
    print_pathway_compact(potato, show_ko = show_ko, show_ec = show_ec)
    return(invisible(potato))
  }

  cli::cli_h1(potato@name)
  cli::cli_text("Source: {.field {potato@source}}")
  cli::cli_text("Tags: {.val {paste(potato@tags, collapse = ', ')}}")
  cli::cli_rule()

  # Group nodes by step
  steps <- list()
  for (node in potato@nodes) {
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
      # Single gene - check if complex
      node <- step_nodes[[1]]
      if (!is.null(node$databases$kofam) && length(node$databases$kofam) > 1) {
        # Complex (multiple KOs in single node)
        cli::cli_text("{.strong COMPLEX} (all subunits required):")
        format_gene(node, show_ko, show_ec)
      } else {
        # Simple single gene
        format_gene(node, show_ko, show_ec)
      }
    } else {
      # Multiple nodes = OR alternatives
      cli::cli_text("{.strong ALTERNATIVES} (any one):")
      for (node in step_nodes) {
        # Check if this individual node is a complex
        if (!is.null(node$databases$kofam) && length(node$databases$kofam) > 1) {
          cli::cli_text("  {.strong COMPLEX}:")
          format_gene(node, show_ko, show_ec, indent = TRUE)
        } else {
          format_gene(node, show_ko, show_ec, indent = TRUE)
        }
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
print_multi_pathway_network <- function(potato, show_ko = FALSE, show_ec = FALSE) {

  cli::cli_h1(potato@name)
  cli::cli_text("Source: {.field {potato@source}}")
  cli::cli_text("Tags: {.val {paste(potato@tags, collapse = ', ')}}")
  cli::cli_text("")
  cli::cli_alert_info("Multi-pathway network with {length(potato@edges)} pathway{?s}")
  cli::cli_rule()

  # Print each pathway
  pathways <- potato@edges  # In network potatoes, edges slot contains pathways

  for (pathway_id in names(pathways)) {
    pathway <- pathways[[pathway_id]]

    cli::cli_h2(pathway$name %||% pathway_id)
    cli::cli_text("Type: {.field {pathway$type}}")
    if (!is.null(pathway$kegg_module)) {
      cli::cli_text("KEGG: {.field {pathway$kegg_module}}")
    }
    if (!is.null(pathway$notes)) {
      cli::cli_text("{.emph {pathway$notes}}")
    }
    cli::cli_text("")

    # Build compact pathway string for this pathway
    pathway_nodes <- pathway$nodes

    # Get steps from pathway-specific nodes
    steps <- list()
    for (node_id in names(pathway_nodes)) {
      node_attrs <- pathway_nodes[[node_id]]
      step_num <- node_attrs$step
      if (is.list(step_num)) step_num <- step_num[[1]]

      if (is.null(steps[[as.character(step_num)]])) {
        steps[[as.character(step_num)]] <- list()
      }

      # Merge global node data with pathway-specific attributes
      global_node <- Find(function(n) n$id == node_id, potato@nodes)
      merged_node <- global_node
      merged_node$step <- node_attrs$step
      merged_node$required <- node_attrs$required
      merged_node$marker <- node_attrs$marker

      steps[[as.character(step_num)]][[length(steps[[as.character(step_num)]]) + 1]] <- merged_node
    }

    # Sort by step
    step_numbers <- as.integer(names(steps))
    steps <- steps[order(step_numbers)]

    # Build compact string
    parts <- character()

    if (!is.null(pathway$input)) {
      parts <- c(parts, paste0("[", pathway$input$compound, "]"))
    }

    for (step_name in names(steps)) {
      step_nodes <- steps[[step_name]]

      if (length(step_nodes) == 1) {
        node <- step_nodes[[1]]
        gene_str <- format_gene_compact(node, show_ko, show_ec)
        parts <- c(parts, gene_str)
      } else {
        # Multiple alternatives
        alt_genes <- sapply(step_nodes, function(node) format_gene_compact(node, show_ko, show_ec))
        parts <- c(parts, paste0("(", paste(alt_genes, collapse = " | "), ")"))
      }
    }

    if (!is.null(pathway$output)) {
      parts <- c(parts, paste0("[", pathway$output$compound, "]"))
    }

    # Print pathway flow
    pathway_str <- paste(parts, collapse = " -> ")
    cli::cli_text("  {pathway_str}")
    cli::cli_text("")
  }

  cli::cli_rule()
  cli::cli_text("{.emph * = marker gene}")
  cli::cli_text("{.emph ^ = optional}")
  cli::cli_text("{.emph {{n}} = complex with n subunits}")
  cli::cli_text("")

  invisible(potato)
}


#' Print compact pathway view (internal)
#' @noRd
print_pathway_compact <- function(potato, show_ko = FALSE, show_ec = TRUE) {

  cli::cli_h2(potato@name)
  cli::cli_text("Source: {.field {potato@source}}")

  # Group nodes by step
  steps <- list()
  for (node in potato@nodes) {
    step_num <- node$step
    if (is.list(step_num)) step_num <- step_num[[1]]

    if (is.null(steps[[as.character(step_num)]])) {
      steps[[as.character(step_num)]] <- list()
    }
    steps[[as.character(step_num)]][[length(steps[[as.character(step_num)]]) + 1]] <- node
  }

  # Sort by step number
  step_numbers <- as.integer(names(steps))
  steps <- steps[order(step_numbers)]

  # Build compact string
  parts <- character()

  # Add input (if field exists and is populated)
  input_val <- tryCatch(potato@input, error = function(e) NULL)
  if (!is.null(input_val) && length(input_val) > 0) {
    parts <- c(parts, paste0("[", input_val$compound, "]"))
  }

  for (step_name in names(steps)) {
    step_nodes <- steps[[step_name]]

    if (length(step_nodes) == 1) {
      # Single gene or complex
      node <- step_nodes[[1]]
      gene_str <- format_gene_compact(node, show_ko, show_ec)
      parts <- c(parts, gene_str)

    } else {
      # Multiple alternatives
      alt_genes <- sapply(step_nodes, function(node) format_gene_compact(node, show_ko, show_ec))
      parts <- c(parts, paste0("(", paste(alt_genes, collapse = " | "), ")"))
    }
  }

  # Add output (if field exists and is populated)
  output_val <- tryCatch(potato@output, error = function(e) NULL)
  if (!is.null(output_val) && length(output_val) > 0) {
    parts <- c(parts, paste0("[", output_val$compound, "]"))
  }

  # Print compact view
  pathway_str <- paste(parts, collapse = " -> ")
  cli::cli_text(pathway_str)
  cli::cli_text("")
  cli::cli_text("{.emph * = marker gene}")
  cli::cli_text("{.emph ^ = optional step}")
  cli::cli_text("{.emph {{n}} = complex with n subunits}")
  if (show_ko) {
    cli::cli_text("{.emph [K#####] = KO identifier}")
  }
  cli::cli_text("")
}


#' Format gene for compact view (internal)
#' @noRd
format_gene_compact <- function(node, show_ko = FALSE, show_ec = TRUE) {
  gene_str <- node$id

  # Add EC if requested
  if (show_ec && !is.null(node$ec) && length(node$ec) > 0) {
    ec_str <- if (length(node$ec) > 1) {
      paste(node$ec, collapse = ",")
    } else {
      node$ec[1]
    }
    gene_str <- paste0(gene_str, "[", ec_str, "]")
  }

  # Add KO if requested (after EC)
  if (show_ko && !is.null(node$databases$kofam)) {
    kos <- unlist(node$databases$kofam)
    if (length(kos) == 1) {
      gene_str <- paste0(gene_str, "[", kos, "]")
    } else if (length(kos) > 1) {
      # Complex - show all KOs with +
      gene_str <- paste0(gene_str, "[", paste(kos, collapse = "+"), "]")
    }
  } else if (!show_ec) {
    # Just show {n} for complex without EC/KO details
    if (!is.null(node$databases$kofam) && length(node$databases$kofam) > 1) {
      gene_str <- paste0(gene_str, "{", length(node$databases$kofam), "}")
    }
  }

  # Add marker indicator
  if (node$marker %||% FALSE) {
    gene_str <- paste0(gene_str, "*")
  }

  # Add optional indicator
  if (!(node$required %||% TRUE)) {
    gene_str <- paste0(gene_str, "^")
  }

  gene_str
}


#' Format a gene for printing (internal)
#' @noRd
format_gene <- function(node, show_ko, show_ec, indent = FALSE) {
  prefix <- if (indent) "  • " else ""

  # Build label
  label <- paste0(prefix, "{.strong ", node$id, "}")

  # Add EC
  if (show_ec && !is.null(node$ec) && length(node$ec) > 0) {
    label <- paste0(label, " {.field [", paste(node$ec, collapse = ", "), "]}")
  }

  # Add name
  label <- paste0(label, " - {.emph ", node$name, "}")

  cli::cli_text(label)

  # Show databases
  if (!is.null(node$databases)) {
    db_info <- character()
    if (!is.null(node$databases$kofam)) {
      kos <- unlist(node$databases$kofam)
      if (show_ko || length(kos) > 1) {  # Always show if complex
        db_info <- c(db_info, paste0("KOfam: ", paste(kos, collapse = "+")))
      }
    }
    if (!is.null(node$databases$blast)) {
      n_blast <- length(unlist(node$databases$blast))
      db_info <- c(db_info, paste0("BLAST: ", n_blast, " refs"))
    }
    if (!is.null(node$databases$hmm)) {
      n_hmm <- length(unlist(node$databases$hmm))
      db_info <- c(db_info, paste0("HMM: ", n_hmm, " profiles"))
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
