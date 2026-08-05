#' Print detailed pathway information
#'
#' Displays detailed information about a specific pathway in a potato,
#' including a markdown table of genes with their properties.
#'
#' @param potato Potato object or path to potato JSON file
#' @param pathway Pathway ID to display details for
#'
#' @return Invisible NULL (prints to console)
#' @export

print_pathway_detail <- function(potato, pathway = NULL) {

  # Load potato if character path
  if (is.character(potato)) {
    potato <- load_potato(potato)
  }

  # Check if multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  if (!is_network) {
    cli::cli_abort("This function only works with multi-pathway network potatoes")
  }

  # Get available pathways with names
  pathway_info <- sapply(names(potato@edges), function(id) {
    name <- potato@edges[[id]]$name %||% id
    paste0(id, " (", name, ")")
  })

  # Check if pathway argument provided
  if (is.null(pathway)) {
    cli::cli_abort(c(
      "Pathway ID required",
      "i" = "Available pathways:",
      paste0("  - ", pathway_info)
    ))
  }

  # Check pathway exists
  if (!pathway %in% names(potato@edges)) {
    cli::cli_abort(c(
      "Pathway {.val {pathway}} not found",
      "i" = "Available pathways:",
      paste0("  - ", pathway_info)
    ))
  }

  pathway_obj <- potato@edges[[pathway]]

  # Print header
  cli::cli_h1(pathway_obj$name %||% pathway)

  # Print metadata
  if (!is.null(pathway_obj$type)) {
    cli::cli_text("{.strong Type:} {pathway_obj$type}")
  }
  if (!is.null(pathway_obj$verified)) {
    verified_status <- if (isTRUE(pathway_obj$verified)) "Yes" else "No"
    cli::cli_text("{.strong Verified:} {verified_status}")
  }
  if (!is.null(pathway_obj$source)) {
    cli::cli_text("{.strong Source:} {pathway_obj$source}")
  }
  if (!is.null(pathway_obj$kegg_module)) {
    cli::cli_text("{.strong KEGG Module:} {pathway_obj$kegg_module}")
  }
  if (!is.null(pathway_obj$notes)) {
    cli::cli_text("")
    cli::cli_text("{pathway_obj$notes}")
  }

  # Print input/output
  cli::cli_text("")
  if (!is.null(pathway_obj$input)) {
    input_str <- pathway_obj$input$compound
    if (!is.null(pathway_obj$input$kegg_compound)) {
      input_str <- paste0(input_str, " (", pathway_obj$input$kegg_compound, ")")
    }
    cli::cli_text("{.strong Input:} {input_str}")
  }
  if (!is.null(pathway_obj$output)) {
    output_str <- pathway_obj$output$compound
    if (!is.null(pathway_obj$output$kegg_compound)) {
      output_str <- paste0(output_str, " (", pathway_obj$output$kegg_compound, ")")
    }
    cli::cli_text("{.strong Output:} {output_str}")
  }

  # Build gene table
  cli::cli_text("")
  cli::cli_h2("Genes")

  # Get pathway-specific node info
  pathway_nodes <- pathway_obj$nodes

  # Get global node info (detection methods)
  gene_rows <- list()

  for (gene_id in names(pathway_nodes)) {
    pathway_node <- pathway_nodes[[gene_id]]

    # Find global node
    global_node <- NULL
    for (node in potato@nodes) {
      if (node$id == gene_id) {
        global_node <- node
        break
      }
    }

    if (is.null(global_node)) {
      next
    }

    # Extract info
    step <- pathway_node$step
    if (is.list(step)) step <- paste(unlist(step), collapse = ", ")

    name <- global_node$name %||% ""

    ec <- if (!is.null(global_node$ec) && length(global_node$ec) > 0) {
      paste(global_node$ec, collapse = ", ")
    } else {
      ""
    }

    # Detection methods
    detection <- character()
    if (!is.null(global_node$databases)) {
      for (db_type in names(global_node$databases)) {
        terms <- global_node$databases[[db_type]]
        if (length(terms) > 0) {
          detection <- c(detection, paste0(db_type, ": ", paste(terms, collapse = ", ")))
        }
      }
    }
    detection_str <- if (length(detection) > 0) paste(detection, collapse = "; ") else ""

    required <- if (isTRUE(pathway_node$required)) "Yes" else "No"
    marker <- if (isTRUE(pathway_node$marker)) "Yes" else "No"

    notes <- global_node$notes %||% ""

    gene_rows[[length(gene_rows) + 1]] <- data.frame(
      Step = step,
      Gene = gene_id,
      Name = name,
      EC = ec,
      Detection = detection_str,
      Required = required,
      Marker = marker,
      Notes = notes,
      stringsAsFactors = FALSE
    )
  }

  # Combine into data frame
  if (length(gene_rows) > 0) {
    gene_table <- do.call(rbind, gene_rows)

    # Sort by step
    gene_table <- gene_table[order(as.numeric(gene_table$Step)), ]

    # Print as markdown table
    cat("\n")
    print(knitr::kable(gene_table, format = "markdown", row.names = FALSE))
    cat("\n")
  }

  # Print scoring info
  cli::cli_text("")
  cli::cli_h2("Scoring")

  if (!is.null(pathway_obj$scoring)) {
    scoring <- pathway_obj$scoring

    if (!is.null(scoring$min_fraction)) {
      cli::cli_text("{.strong Threshold:} {scoring$min_fraction} ({scoring$min_fraction * 100}% of steps required)")
    }
    if (!is.null(scoring$marker_mode)) {
      cli::cli_text("{.strong Marker mode:} {scoring$marker_mode}")
    }
    if (!is.null(scoring$notes)) {
      cli::cli_text("")
      cli::cli_text("{scoring$notes}")
    }
  }

  invisible(NULL)
}
