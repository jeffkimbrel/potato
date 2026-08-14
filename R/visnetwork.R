#' Update potato JSON with node coordinates
#'
#' After arranging nodes in visNetwork and exporting coordinates,
#' use this function to add x,y fields to your potato JSON.
#'
#' @param potato_path Path to potato JSON file
#' @param coords_path Path to coordinates JSON file (from visNetwork export)
#' @param output_path Path to save updated potato JSON (default: overwrites original)
#' @param with_compounds Logical - if NULL (default), auto-detects from coordinate file; if TRUE/FALSE, forces that mode
#'
#' @export
update_potato_coordinates <- function(potato_path, coords_path, output_path = NULL, with_compounds = NULL) {

  if (is.null(output_path)) {
    output_path <- potato_path
  }

  # Read and validate potato JSON
  potato_json <- jsonlite::read_json(potato_path, simplifyVector = FALSE)

  # Validate potato structure
  validation <- validate_potato(potato_json, strict = FALSE)
  if (!validation$valid) {
    cli::cli_abort(c(
      "Potato validation failed: {basename(potato_path)}",
      "x" = "Errors found:",
      set_names(validation$errors, rep("*", length(validation$errors)))
    ))
  }

  # Read coordinates
  coords <- jsonlite::read_json(coords_path, simplifyVector = TRUE)

  # Normalize coordinates to reasonable range for plotting
  # visNetwork uses pixel coordinates (e.g., -500 to 500)
  # Static plots work better with normalized coords (e.g., -20 to 20)
  x_range <- range(coords$x)
  y_range <- range(coords$y)
  x_span <- diff(x_range)
  y_span <- diff(y_range)

  # Scale to roughly -20 to 20 range, maintaining aspect ratio
  max_span <- max(x_span, y_span)
  scale_factor <- 40 / max_span  # Target range of 40 units (-20 to 20)

  coords$x <- (coords$x - mean(x_range)) * scale_factor
  coords$y <- (coords$y - mean(y_range)) * scale_factor

  cli::cli_alert_info("Normalized coordinates: X [{round(min(coords$x), 2)}, {round(max(coords$x), 2)}], Y [{round(min(coords$y), 2)}, {round(max(coords$y), 2)}]")

  # Separate enzyme nodes and compound nodes (including INPUT/OUTPUT)
  enzyme_nodes <- coords[!grepl("^(COMPOUND_|INPUT_|OUTPUT_)", coords$id), ]
  compound_nodes <- coords[grepl("^(COMPOUND_|INPUT_|OUTPUT_)", coords$id), ]

  # Auto-detect if coordinates include compounds (if not explicitly specified)
  if (is.null(with_compounds)) {
    has_compounds <- nrow(compound_nodes) > 0
    with_compounds <- has_compounds
    if (has_compounds) {
      cli::cli_alert_info("Auto-detected: coordinates include compound nodes ({nrow(compound_nodes)} compound nodes found)")
    } else {
      cli::cli_alert_info("Auto-detected: coordinates are enzyme-only ({nrow(enzyme_nodes)} enzyme nodes)")
    }
  }

  # Show which genes will be updated
  genes_in_potato <- sapply(potato_json$genes, function(g) g$id)
  genes_in_coords <- enzyme_nodes$id
  missing_in_coords <- setdiff(genes_in_potato, genes_in_coords)
  extra_in_coords <- setdiff(genes_in_coords, genes_in_potato)

  if (length(missing_in_coords) > 0) {
    cli::cli_alert_warning("Genes in potato but not in coordinates: {paste(missing_in_coords, collapse=', ')}")
  }
  if (length(extra_in_coords) > 0) {
    cli::cli_alert_warning("Genes in coordinates but not in potato: {paste(extra_in_coords, collapse=', ')}")
  }

  # Create lookup table for gene coordinates
  gene_coords_lookup <- stats::setNames(
    split(enzyme_nodes[, c("x", "y")], seq_len(nrow(enzyme_nodes))),
    enzyme_nodes$id
  )

  # Determine field names based on with_compounds
  x_field <- if (with_compounds) "x_compounds" else "x"
  y_field <- if (with_compounds) "y_compounds" else "y"

  # Update gene coordinates
  for (i in seq_along(potato_json$genes)) {
    gene_id <- potato_json$genes[[i]]$id

    if (gene_id %in% names(gene_coords_lookup)) {
      potato_json$genes[[i]][[x_field]] <- gene_coords_lookup[[gene_id]]$x
      potato_json$genes[[i]][[y_field]] <- gene_coords_lookup[[gene_id]]$y
    }
  }

  # Update compound node coordinates (stored separately in compound_coordinates)
  if (nrow(compound_nodes) > 0) {
    # Create/update compound_coordinates field
    if (is.null(potato_json$compound_coordinates)) {
      potato_json$compound_coordinates <- list()
    }

    # Store compound coordinates (includes COMPOUND_, INPUT_, OUTPUT_)
    for (i in seq_len(nrow(compound_nodes))) {
      compound_id <- compound_nodes$id[i]
      potato_json$compound_coordinates[[compound_id]] <- list(
        id = compound_id,
        x = compound_nodes$x[i],
        y = compound_nodes$y[i]
      )
    }

    cli::cli_alert_info("Saved {nrow(compound_nodes)} compound node coordinates (includes INPUT/OUTPUT)")
  }

  # Write updated JSON
  jsonlite::write_json(
    potato_json,
    output_path,
    pretty = TRUE,
    auto_unbox = TRUE
  )

  coord_type <- if (with_compounds) "with compounds" else "without compounds"
  cli::cli_alert_success("Updated potato JSON with coordinates ({coord_type}): {.file {output_path}}")
  cli::cli_alert_info("Total nodes updated: {sum(sapply(potato_json$nodes, function(n) !is.null(n[[x_field]])))} using fields {.field {x_field}}/{.field {y_field}}")

  # Provide usage guidance
  if (with_compounds) {
    cli::cli_alert_info("To use these coordinates, plot with: {.code show_compounds = TRUE}")
  } else {
    cli::cli_alert_info("To use these coordinates, plot with: {.code show_compounds = FALSE} or omit the parameter")
  }

  invisible(output_path)
}
