#' View pathway details as HTML table in browser
#'
#' Opens a simple HTML table in the system browser showing pathway gene details.
#' Useful when terminal output is too wide to read comfortably.
#'
#' @param potato Potato S7 object or path to JSON
#' @param pathway For multi-pathway networks, which pathway to show (single pathway ID only)
#' @param layout Layout algorithm for visualization: "xy" (curated coords), "fr", "kk", "sugiyama", "tree", "circle", "grid" (default: "fr")
#'
#' @export
view_pathway_detail <- function(potato, pathway = NULL, layout = "fr") {

  # Load if path provided (V2 only)
  if (is.character(potato)) {
    potato <- load_potato_v2(potato)  # Validates automatically
  }

  # Check if multi-pathway network (> 1 pathway)
  is_multi_pathway <- !is.null(potato@pathways) && length(potato@pathways) > 1

  # Get available pathways
  available_pathways <- names(potato@pathways)

  # For single-pathway, set pathway to first one
  if (is.null(pathway)) {
    if (is_multi_pathway) {
      cli::cli_abort(c(
        "Multi-pathway network requires pathway parameter",
        "i" = "Available pathways: {paste(available_pathways, collapse=', ')}"
      ))
    }
    pathway <- names(potato@pathways)[1]
  }

  if (!is.null(pathway) && length(pathway) > 1) {
    cli::cli_abort(c(
      "view_pathway_detail() only supports a single pathway",
      "i" = "For plotting multiple pathways, use plot_potato_static() or plot_potato_interactive() with pathway = c(...)"
    ))
  }

  if (!is.null(pathway) && !pathway %in% available_pathways) {
    available <- paste(available_pathways, collapse = ", ")
    cli::cli_abort("Pathway '{pathway}' not found. Available: {available}")
  }

  # Build rows (V2 schema)
  pathway_info <- potato@pathways[[pathway]]
  title <- paste0(potato@name, " - ", pathway_info$name %||% pathway)

  # Get genes used in this pathway
  # Option 1: Explicit genes array in pathway (includes non-catalytic genes)
  # Option 2: Extract from edges (catalytic genes only)
  # Get all gene IDs from genes array AND edges (to capture complexes)
  pathway_gene_ids <- character()

  if (!is.null(pathway_info$genes) && length(pathway_info$genes) > 0) {
    pathway_gene_ids <- pathway_info$genes
  }

  # Also extract from edges (may include complex IDs not in genes array)
  edge_ids <- unique(unlist(lapply(pathway_info$edges, function(e) {
    c(e$from, e$to)
  })))
  # Filter to only genes (exclude compound IDs)
  all_gene_ids <- sapply(potato@genes, function(g) g$id)
  edge_gene_ids <- edge_ids[edge_ids %in% all_gene_ids]

  # Combine: components from genes array + complexes from edges
  pathway_gene_ids <- unique(c(pathway_gene_ids, edge_gene_ids))

  # Filter genes to only those in this pathway (includes components AND complexes)
  pathway_genes <- Filter(function(g) g$id %in% pathway_gene_ids, potato@genes)

  # Helper to find which pathways contain a gene
  get_pathways_for_gene <- function(gene_id) {
    pathway_names <- c()
    for (pw_id in names(potato@pathways)) {
      pw <- potato@pathways[[pw_id]]
      # Check if gene appears in any edge
      gene_in_pathway <- any(sapply(pw$edges, function(e) {
        e$from == gene_id || e$to == gene_id
      }))
      if (gene_in_pathway) {
        pathway_names <- c(pathway_names, pw_id)
      }
    }
    pathway_names
  }

  rows <- lapply(pathway_genes, function(gene) {
    # Check if this is a complex
    is_complex <- !is.null(gene$type) && gene$type == "complex"

    if (is_complex) {
      # Complex entry - show components instead of databases
      component_str <- paste(gene$components, collapse = ", ")
      db_str <- paste0("<em>Complex components: ", component_str, "</em>")

      # Get pathway-specific attributes from edges
      pathway_edges <- pathway_info$edges
      is_required <- any(sapply(pathway_edges, function(e) {
        (e$from == gene$id || e$to == gene$id) && !is.null(e$required) && e$required
      }))
      is_marker <- any(sapply(pathway_edges, function(e) {
        (e$from == gene$id || e$to == gene$id) && !is.null(e$marker) && e$marker
      }))

      # Get all pathways this complex appears in
      gene_pathways <- get_pathways_for_gene(gene$id)
      pathway_display <- sapply(gene_pathways, function(pw) {
        if (pw == pathway) {
          paste0("<strong>", pw, "</strong>")
        } else {
          pw
        }
      })

      list(
        id = gene$id,
        required = is_required,
        marker = is_marker,
        name = gene$name,
        ec = "",
        reactions = "",
        databases = db_str,
        pathways = paste(pathway_display, collapse = ", "),
        notes = gene$notes %||% ""
      )
    } else {
      # Regular gene entry
      # Format databases
      dbs <- sapply(names(gene$databases), function(db) {
        terms <- paste(gene$databases[[db]], collapse = ", ")
        paste0(db, ": ", terms)
      })
      db_str <- paste(dbs, collapse = "<br>")

      # Get pathway-specific attributes from edges
      pathway_edges <- pathway_info$edges
      is_required <- any(sapply(pathway_edges, function(e) {
        (e$from == gene$id || e$to == gene$id) && !is.null(e$required) && e$required
      }))
      is_marker <- any(sapply(pathway_edges, function(e) {
        (e$from == gene$id || e$to == gene$id) && !is.null(e$marker) && e$marker
      }))

      # Get reactions used in pathway for this gene
      pathway_reactions <- unique(sapply(pathway_edges, function(e) {
        if ((e$from == gene$id || e$to == gene$id) && !is.null(e$reaction)) {
          e$reaction
        } else {
          NA
        }
      }))
      pathway_reactions <- pathway_reactions[!is.na(pathway_reactions)]

      # Format reactions with color coding
      reactions_formatted <- if (length(gene$reactions) > 0) {
        sapply(gene$reactions, function(rxn) {
          if (rxn %in% pathway_reactions) {
            # Pathway reaction - normal color
            rxn
          } else {
            # Not in pathway - gray
            paste0("<span style='color: #999;'>", rxn, "</span>")
          }
        })
      } else {
        character(0)
      }

      # Get all pathways this gene appears in
      gene_pathways <- get_pathways_for_gene(gene$id)
      # Highlight current pathway in bold
      pathway_display <- sapply(gene_pathways, function(pw) {
        if (pw == pathway) {
          paste0("<strong>", pw, "</strong>")
        } else {
          pw
        }
      })

      list(
        id = gene$id,
        required = is_required,
        marker = is_marker,
        name = gene$name,
        ec = paste(gene$ec, collapse = ", "),
        reactions = paste(reactions_formatted, collapse = ", "),
        databases = db_str,
        pathways = paste(pathway_display, collapse = ", "),
        notes = gene$notes %||% ""
      )
    }
  })

  # Generate plot
  plot_widget <- NULL
  tryCatch({
    # Create interactive plot with fixed height for embedding
    g <- build_graph_v2(potato, pathway_id = pathway)
    plot_widget <- plot_v2_interactive(g, layout = layout, height = "600px")
  }, error = function(e) {
    cli::cli_warn("Could not generate plot: {e$message}")
  })

  # Get pathway-specific metadata (V2 schema)
  source_info <- potato@source
  potato_notes <- potato@notes %||% ""
  pathway_notes <- pathway_info$notes %||% ""

  # Extract input/output info (array of compound IDs)
  input_info <- if (!is.null(pathway_info$input)) {
    compound_names <- sapply(pathway_info$input, function(cid) {
      comp <- Filter(function(c) c$id == cid, potato@compounds)
      if (length(comp) > 0) comp[[1]]$name else cid
    })
    list(display = paste(compound_names, collapse = ", "))
  } else NULL

  output_info <- if (!is.null(pathway_info$output)) {
    compound_names <- sapply(pathway_info$output, function(cid) {
      comp <- Filter(function(c) c$id == cid, potato@compounds)
      if (length(comp) > 0) comp[[1]]$name else cid
    })
    list(display = paste(compound_names, collapse = ", "))
  } else NULL

  scoring_info <- pathway_info$scoring
  verified <- !is.null(pathway_info$verified) && pathway_info$verified == TRUE

  # Build verification status banner
  verification_html <- if (verified) {
    "<div style='background: #d4edda; padding: 15px; margin-bottom: 15px; border-left: 4px solid #28a745; color: #155724;'>
    <strong>✓ VERIFIED</strong> - This pathway has been manually validated.
    </div>"
  } else {
    "<div style='background: #fff3cd; padding: 15px; margin-bottom: 15px; border-left: 4px solid #ffc107; color: #856404;'>
    <strong>⚠ UNVERIFIED</strong> - This pathway has not been validated.
    </div>"
  }

  # Build metadata section
  metadata_html <- paste0(
    verification_html,
    "<div style='background: #f9f9f9; padding: 15px; margin-bottom: 20px; border-left: 4px solid #4CAF50;'>",
    "<p><strong>Source:</strong> ", source_info, "</p>",
    if (nchar(potato_notes) > 0) paste0("<p><strong>Potato Description:</strong> ", potato_notes, "</p>") else "",
    if (nchar(pathway_notes) > 0) paste0("<p><strong>Pathway Notes:</strong> ", pathway_notes, "</p>") else "",
    if (!is.null(input_info)) paste0("<p><strong>Input:</strong> ", input_info$display, "</p>") else "",
    if (!is.null(output_info)) paste0("<p><strong>Output:</strong> ", output_info$display, "</p>") else "",
    "</div>"
  )

  # Build gene table
  table_rows <- sapply(rows, function(row) {
    paste0(
      "<tr>",
      "<td>", row$id, "</td>",
      "<td>", ifelse(row$required, "✓", ""), "</td>",
      "<td>", ifelse(row$marker, "✓", ""), "</td>",
      "<td>", row$name, "</td>",
      "<td>", row$ec, "</td>",
      "<td>", row$reactions, "</td>",
      "<td>", row$databases, "</td>",
      "<td>", row$pathways, "</td>",
      "<td style='max-width: 400px;'>", row$notes, "</td>",
      "</tr>"
    )
  })

  # Build edges section (V2 schema)
  edges_data <- pathway_info$edges

  # Helper to look up compound name from ID, or format gene with blue color
  get_node_label <- function(node_id) {
    # V2: search compounds list
    compound <- Find(function(c) c$id == node_id, potato@compounds)
    if (!is.null(compound)) {
      # It's a compound - return formatted name (no color)
      return(paste0(compound$name, " [", node_id, "]"))
    }
    # It's a gene - return with blue color
    return(paste0("<span style='color: #216ADE;'>", node_id, "</span>"))
  }

  # Build consolidated rows: compound → gene → compound
  edges_rows <- if (length(edges_data) > 0 && length(pathway_gene_ids) > 0) {

    # Get compound IDs for checking
    compound_ids <- sapply(potato@compounds, function(c) c$id)

    # Order genes by when their input compound first appears in the pathway
    # Walk through edges and assign order based on first appearance
    gene_order <- setNames(rep(Inf, length(pathway_gene_ids)), pathway_gene_ids)
    current_order <- 1

    # Track which compounds have been produced
    compounds_seen <- character(0)

    # Add input compounds as already seen
    if (!is.null(pathway_info$input)) {
      compounds_seen <- c(compounds_seen, pathway_info$input)
    }

    # Iterate until all genes are ordered
    max_iterations <- length(pathway_gene_ids) * 2
    iteration <- 0

    while (any(is.infinite(gene_order)) && iteration < max_iterations) {
      iteration <- iteration + 1

      for (gene_id in pathway_gene_ids) {
        if (!is.infinite(gene_order[gene_id])) next  # Already ordered

        # Find input compounds for this gene
        input_compounds <- sapply(edges_data, function(e) {
          if (e$to == gene_id && e$from %in% compound_ids) e$from else NA
        })
        input_compounds <- input_compounds[!is.na(input_compounds)]

        # If any input compound is available, order this gene
        if (length(input_compounds) > 0 && any(input_compounds %in% compounds_seen)) {
          gene_order[gene_id] <- current_order
          current_order <- current_order + 1

          # Mark this gene's output compounds as seen
          output_compounds <- sapply(edges_data, function(e) {
            if (e$from == gene_id && e$to %in% compound_ids) e$to else NA
          })
          output_compounds <- output_compounds[!is.na(output_compounds)]
          compounds_seen <- c(compounds_seen, output_compounds)
        }
      }
    }

    # Sort genes by order
    sorted_gene_ids <- names(sort(gene_order))

    # For each gene, find its input and output compounds
    unlist(lapply(sorted_gene_ids, function(gene_id) {
      # Find edges going into this gene (compound → gene)
      input_edges <- Filter(function(e) e$to == gene_id, edges_data)
      # Find edges going out of this gene (gene → compound)
      output_edges <- Filter(function(e) e$from == gene_id, edges_data)

      # Create a row for each input-gene-output combination
      if (length(input_edges) > 0 && length(output_edges) > 0) {
        unlist(lapply(input_edges, function(in_edge) {
          lapply(output_edges, function(out_edge) {
            input_label <- get_node_label(in_edge$from)
            gene_label <- get_node_label(gene_id)
            output_label <- get_node_label(out_edge$to)

            # Add reaction if present (use output edge reaction)
            reaction_str <- if (!is.null(out_edge$reaction)) {
              paste0(" <strong>[", out_edge$reaction, "]</strong>")
            } else ""

            paste0("<li>", input_label, " → ", gene_label, " → ", output_label, reaction_str, "</li>")
          })
        }))
      } else if (length(input_edges) > 0) {
        # Gene has inputs but no outputs (terminal gene)
        lapply(input_edges, function(in_edge) {
          input_label <- get_node_label(in_edge$from)
          gene_label <- get_node_label(gene_id)
          paste0("<li>", input_label, " → ", gene_label, "</li>")
        })
      } else if (length(output_edges) > 0) {
        # Gene has outputs but no inputs (initial gene)
        lapply(output_edges, function(out_edge) {
          gene_label <- get_node_label(gene_id)
          output_label <- get_node_label(out_edge$to)
          paste0("<li>", gene_label, " → ", output_label, "</li>")
        })
      } else {
        NULL
      }
    }))
  } else {
    "<li><em>No edges defined</em></li>"
  }

  # Handle case where no rows were generated
  if (is.null(edges_rows) || length(edges_rows) == 0) {
    edges_rows <- "<li><em>No edges defined</em></li>"
  }

  # Remove duplicates
  edges_rows <- unique(edges_rows)

  edges_html <- paste0(
    "<h3>Pathway Topology</h3>",
    "<ul style='list-style-type: none; padding-left: 0;'>",
    paste(edges_rows, collapse = "\n"),
    "</ul>"
  )

  # Build scoring section
  scoring_html <- paste0(
    "<h3>Scoring Parameters</h3>",
    "<div style='background: #f0f0f0; padding: 10px; margin-bottom: 20px;'>",
    if (!is.null(scoring_info$min_fraction)) paste0("<p><strong>Min Fraction:</strong> ", scoring_info$min_fraction, "</p>") else "",
    if (!is.null(scoring_info$max_gaps)) paste0("<p><strong>Max Gaps:</strong> ", scoring_info$max_gaps, "</p>") else "",
    "<p><strong>Marker Mode:</strong> ", scoring_info$marker_mode %||% "any", "</p>",
    "</div>"
  )


  # If we have a widget, save it separately and build the full page
  if (!is.null(plot_widget)) {
    # Build HTML page with widget embedded
    page <- htmltools::tagList(
      htmltools::tags$head(
        htmltools::tags$meta(charset = "utf-8"),
        htmltools::tags$title(title),
        htmltools::tags$style(htmltools::HTML("
          body { font-family: sans-serif; margin: 20px; max-width: 1600px; }
          table { border-collapse: collapse; width: 100%; margin-bottom: 20px; table-layout: auto; }
          th, td { border: 1px solid #ddd; padding: 8px; text-align: left; vertical-align: top; word-wrap: break-word; }
          th { background-color: #f2f2f2; font-weight: bold; }
          tr:hover { background-color: #f5f5f5; }
          h2 { color: #333; }
          h3 { color: #666; margin-top: 20px; }
        "))
      ),
      htmltools::tags$body(
        htmltools::tags$h1(title),
        htmltools::HTML(metadata_html),
        htmltools::tags$h3("Pathway Visualization"),
        htmltools::tags$div(style = "margin: 20px 0; height: 600px;", plot_widget),
        htmltools::tags$h3("Genes"),
        htmltools::HTML(paste0(
          "<table><tr>",
          "<th>Gene ID</th><th>Required</th><th>Marker</th>",
          "<th>Name</th><th>EC</th><th>Reactions</th><th>Databases</th><th>Pathways</th><th>Notes</th>",
          "</tr>",
          paste(table_rows, collapse = "\n"),
          "</table>"
        )),
        htmltools::HTML(edges_html),
        htmltools::HTML(scoring_html)
      )
    )

    temp_file <- tempfile(fileext = ".html")
    htmltools::save_html(page, temp_file)
  } else {
    # No widget - use simple HTML
    html <- paste0(
      "<!DOCTYPE html>",
      "<html>",
      "<head>",
      "<meta charset='utf-8'>",
      "<title>", title, "</title>",
      "<style>",
      "body { font-family: sans-serif; margin: 20px; max-width: 1600px; }",
      "table { border-collapse: collapse; width: 100%; margin-bottom: 20px; table-layout: auto; }",
      "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; vertical-align: top; word-wrap: break-word; }",
      "th { background-color: #f2f2f2; font-weight: bold; }",
      "tr:hover { background-color: #f5f5f5; }",
      "h2 { color: #333; }",
      "h3 { color: #666; margin-top: 20px; }",
      "</style>",
      "</head>",
      "<body>",
      "<h1>", title, "</h1>",
      metadata_html,
      "<h3>Genes</h3>",
      "<table>",
      "<tr>",
      "<th>Gene ID</th>",
      "<th>Required</th>",
      "<th>Marker</th>",
      "<th>Name</th>",
      "<th>EC</th>",
      "<th>Reactions</th>",
      "<th>Databases</th>",
      "<th>Pathways</th>",
      "<th>Notes</th>",
      "</tr>",
      paste(table_rows, collapse = "\n"),
      "</table>",
      edges_html,
      scoring_html,
      "</body>",
      "</html>"
    )

    temp_file <- tempfile(fileext = ".html")
    writeLines(html, temp_file)
  }

  cli::cli_alert_success("Opening in browser: {.file {basename(temp_file)}}")
  browseURL(temp_file)

  invisible(temp_file)
}
