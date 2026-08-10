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

  # Load if path provided
  if (is.character(potato)) {
    potato <- load_potato(potato)
  }

  # Check if multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  if (is_network && is.null(pathway)) {
    cli::cli_abort(c(
      "Multi-pathway network requires pathway parameter",
      "i" = "Available pathways: {paste(names(potato@edges), collapse=', ')}"
    ))
  }

  if (is_network && !is.null(pathway)) {
    # Multi-pathway: extract specific pathway (single only for table view)
    if (length(pathway) > 1) {
      cli::cli_abort(c(
        "view_pathway_detail() only supports a single pathway",
        "i" = "For plotting multiple pathways, use plot_potato_static() or plot_potato_interactive() with pathway = c(...)"
      ))
    }

    if (!pathway %in% names(potato@edges)) {
      available <- paste(names(potato@edges), collapse = ", ")
      cli::cli_abort("Pathway '{pathway}' not found. Available: {available}")
    }

    pathway_info <- potato@edges[[pathway]]
    pathway_nodes <- pathway_info$nodes
    title <- paste0(potato@name, " - ", pathway_info$name %||% pathway)

    # Build rows
    rows <- lapply(names(pathway_nodes), function(node_id) {
      global_node <- Find(function(n) n$id == node_id, potato@nodes)
      pathway_node <- pathway_nodes[[node_id]]

      # Format databases
      dbs <- sapply(names(global_node$databases), function(db) {
        terms <- paste(global_node$databases[[db]], collapse = ", ")
        paste0(db, ": ", terms)
      })
      db_str <- paste(dbs, collapse = "<br>")

      list(
        id = node_id,
        step = pathway_node$step,
        required = pathway_node$required,
        marker = pathway_node$marker,
        name = global_node$name,
        ec = paste(global_node$ec, collapse = ", "),
        databases = db_str,
        notes = global_node$notes %||% ""
      )
    })

  } else {
    # Single-pathway potato
    title <- potato@name

    rows <- lapply(potato@nodes, function(node) {
      # Format databases
      dbs <- sapply(names(node$databases), function(db) {
        terms <- paste(node$databases[[db]], collapse = ", ")
        paste0(db, ": ", terms)
      })
      db_str <- paste(dbs, collapse = "<br>")

      list(
        id = node$id,
        step = node$step,
        required = node$required %||% TRUE,
        marker = node$marker %||% FALSE,
        name = node$name,
        ec = paste(node$ec, collapse = ", "),
        databases = db_str,
        notes = node$notes %||% ""
      )
    })
  }

  # Generate plot
  plot_img <- NULL
  tryCatch({
    # Create plot
    if (is_network) {
      p <- plot_potato_static(potato, pathway = pathway, show_hulls = FALSE, show_compounds = TRUE, layout = layout)
    } else {
      p <- plot_potato_static(potato, show_hulls = FALSE, show_compounds = TRUE, layout = layout)
    }

    # Save to temp PNG
    temp_plot <- tempfile(fileext = ".png")
    ggplot2::ggsave(temp_plot, plot = p, width = 8, height = 6, dpi = 150, bg = "white")

    # Read and encode as base64
    img_raw <- readBin(temp_plot, "raw", file.info(temp_plot)$size)
    plot_img <- paste0("data:image/png;base64,", base64enc::base64encode(img_raw))

    # Clean up temp plot file
    unlink(temp_plot)
  }, error = function(e) {
    cli::cli_warn("Could not generate plot: {e$message}")
  })

  # Get pathway-specific metadata
  if (is_network) {
    source_info <- paste0("KEGG: ", pathway_info$kegg_module %||% "N/A")
    notes_info <- pathway_info$notes %||% ""
    input_info <- pathway_info$input
    output_info <- pathway_info$output
    scoring_info <- pathway_info$scoring
    verified <- !is.null(pathway_info$verified) && pathway_info$verified == TRUE
  } else {
    source_info <- potato@source
    notes_info <- potato@notes %||% ""
    input_info <- tryCatch(potato@input, error = function(e) NULL)
    output_info <- tryCatch(potato@output, error = function(e) NULL)
    scoring_info <- potato@scoring
    # Check for verified field in top-level data
    verified <- FALSE
    if (!is.null(potato@json_path) && file.exists(potato@json_path)) {
      data <- jsonlite::read_json(potato@json_path, simplifyVector = FALSE)
      verified <- !is.null(data$verified) && data$verified == TRUE
    }
  }

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
      "<td>", row$databases, "</td>",
      "<td style='max-width: 400px;'>", row$notes, "</td>",
      "</tr>"
    )
  })

  # Build edges section
  if (is_network) {
    edges_data <- pathway_info$edges
  } else {
    edges_data <- potato@edges
  }

  edges_rows <- if (length(edges_data) > 0) {
    sapply(edges_data, function(edge) {
      compound_str <- if (!is.null(edge$compound)) {
        paste0(" <em>(", edge$compound, ")</em>")
      } else ""
      paste0("<li>", edge$from, " → ", edge$to, compound_str, "</li>")
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

  # Build plot section if available
  plot_html <- if (!is.null(plot_img)) {
    paste0(
      "<h3>Pathway Visualization</h3>",
      "<div style='text-align: center; margin: 20px 0;'>",
      "<img src='", plot_img, "' style='max-width: 100%; height: auto; border: 1px solid #ddd;'>",
      "</div>"
    )
  } else ""

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
    plot_html,
    "<h3>Genes</h3>",
    "<table>",
    "<tr>",
    "<th>Gene ID</th>",
    "<th>Step</th>",
    "<th>Required</th>",
    "<th>Marker</th>",
    "<th>Name</th>",
    "<th>EC</th>",
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

  # Write to temp file and open
  temp_file <- tempfile(fileext = ".html")
  writeLines(html, temp_file)

  cli::cli_alert_success("Opening in browser: {.file {basename(temp_file)}}")
  browseURL(temp_file)

  invisible(temp_file)
}
