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

  rows <- lapply(potato@genes, function(gene) {
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

    list(
      id = gene$id,
      step = NA,  # V2 doesn't have step numbers
      required = is_required,
      marker = is_marker,
      name = gene$name,
      ec = paste(gene$ec, collapse = ", "),
      reactions = paste(reactions_formatted, collapse = ", "),
      databases = db_str,
      notes = gene$notes %||% ""
    )
  })

  # Generate plot
  plot_widget <- NULL
  tryCatch({
    # Create interactive plot with fixed height for embedding
    g <- build_graph_v2(potato)
    plot_widget <- plot_v2_interactive(g, layout = layout, height = "600px")
  }, error = function(e) {
    cli::cli_warn("Could not generate plot: {e$message}")
  })

  # Get pathway-specific metadata (V2 schema)
  source_info <- potato@source
  notes_info <- pathway_info$notes %||% ""
  input_info <- NULL  # V2 doesn't have input/output at pathway level yet
  output_info <- NULL
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
    if (nchar(notes_info) > 0) paste0("<p><strong>Notes:</strong> ", notes_info, "</p>") else "",
    if (!is.null(input_info)) paste0("<p><strong>Input:</strong> ", input_info$compound, " (", input_info$kegg_compound %||% "", ")</p>") else "",
    if (!is.null(output_info)) paste0("<p><strong>Output:</strong> ", output_info$compound, " (", output_info$kegg_compound %||% "", ")</p>") else "",
    "</div>"
  )

  # Build gene table
  table_rows <- sapply(rows, function(row) {
    paste0(
      "<tr>",
      "<td>", row$id, "</td>",
      "<td>", row$step, "</td>",
      "<td>", ifelse(row$required, "✓", ""), "</td>",
      "<td>", ifelse(row$marker, "✓", ""), "</td>",
      "<td>", row$name, "</td>",
      "<td>", row$ec, "</td>",
      "<td>", row$reactions, "</td>",
      "<td>", row$databases, "</td>",
      "<td style='max-width: 400px;'>", row$notes, "</td>",
      "</tr>"
    )
  })

  # Build edges section (V2 schema)
  edges_data <- pathway_info$edges

  # Helper to look up compound name from ID
  get_compound_name <- function(node_id) {
    # V2: search compounds list
    compound <- Find(function(c) c$id == node_id, potato@compounds)
    if (!is.null(compound)) {
      return(paste0(compound$name, " [", node_id, "]"))
    }
    # Return ID if not found or not a compound
    return(node_id)
  }

  edges_rows <- if (length(edges_data) > 0) {
    sapply(edges_data, function(edge) {
      # Convert IDs to names for compounds
      from_label <- get_compound_name(edge$from)
      to_label <- get_compound_name(edge$to)

      # Add reaction if present
      reaction_str <- if (!is.null(edge$reaction)) {
        paste0(" <strong>[", edge$reaction, "]</strong>")
      } else ""

      paste0("<li>", from_label, " → ", to_label, reaction_str, "</li>")
    })
  } else {
    "<li><em>No edges defined</em></li>"
  }

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
    "<p><strong>Min Fraction:</strong> ", scoring_info$min_fraction %||% 0.75, "</p>",
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
          "<th>Gene ID</th><th>Step</th><th>Required</th><th>Marker</th>",
          "<th>Name</th><th>EC</th><th>Reactions</th><th>Databases</th><th>Notes</th>",
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
      "<th>Step</th>",
      "<th>Required</th>",
      "<th>Marker</th>",
      "<th>Name</th>",
      "<th>EC</th>",
      "<th>Reactions</th>",
      "<th>Databases</th>",
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
