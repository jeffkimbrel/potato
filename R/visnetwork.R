#' Build visNetwork interactive plot from potato (internal)
#' @noRd
build_visnetwork <- function(potato, g, node_coords, node_status, has_genome) {

  # Check if multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  # Check if using curated coordinates (already in screen space, typically 100-1000 range)
  # Layout algorithms produce coordinates in range ~-5 to 5, so we scale those up
  using_curated_coords <- any(!is.na(node_coords$x)) &&
                          max(abs(node_coords$x), na.rm = TRUE) > 50

  # Build nodes data frame for visNetwork
  # Curated coords: use as-is (they're already in screen space)
  # Layout coords: scale by 100 to match screen space
  nodes_df <- data.frame(
    id = node_coords$name,
    label = node_status$label,
    title = node_status$hover_text,  # Tooltip on hover
    x = if (using_curated_coords) node_coords$x else node_coords$x * 100,
    y = if (using_curated_coords) node_coords$y else node_coords$y * 100,
    shape = ifelse(node_status$node_shape == "circle", "dot",
            ifelse(node_status$node_shape == "triangle", "triangle",
            ifelse(node_status$node_shape == "diamond", "diamond", "square"))),
    stringsAsFactors = FALSE
  )

  # Add pathway group membership for multi-pathway networks
  if (is_network) {
    # Assign each node to pathway groups
    node_groups <- list()
    for (node_id in nodes_df$id) {
      # Skip compound nodes
      if (grepl("^COMPOUND_", node_id)) {
        node_groups[[node_id]] <- "Compound"
        next
      }

      # Find which pathways contain this gene
      pathways_containing <- character()
      for (pathway_id in names(potato@edges)) {
        pathway <- potato@edges[[pathway_id]]
        if (node_id %in% names(pathway$nodes)) {
          pathways_containing <- c(pathways_containing, pathway$name %||% pathway_id)
        }
      }

      # If node is in multiple pathways, join with " & "
      if (length(pathways_containing) > 0) {
        node_groups[[node_id]] <- paste(pathways_containing, collapse = " & ")
      } else {
        node_groups[[node_id]] <- "Other"
      }
    }

    nodes_df$group <- unlist(node_groups[nodes_df$id])
  }

  # Add colors based on detection status (fill)
  if (has_genome) {
    nodes_df$color.background <- sapply(seq_len(nrow(nodes_df)), function(i) {
      if (node_status$is_compound_node[i]) {
        # Color by compound type using status
        if (node_status$status[i] == "Input") {
          return("#7FD399")  # Green for inputs
        } else if (node_status$status[i] == "Output") {
          return("#F77370")  # Red for outputs
        } else {
          return("#B3B3B3")  # Gray for intermediates
        }
      } else if (node_status$status[i] == "Detected") {
        return("#4CAF50")
      } else if (node_status$status[i] == "Partial") {
        return("#FFA726")
      } else if (node_status$status[i] == "Not detected") {
        return("#F44336")
      } else {
        return("#2196F3")  # Unknown
      }
    })
  } else {
    # No genome - color compounds by type (input=green, output=red, intermediate=gray)
    nodes_df$color.background <- sapply(seq_len(nrow(nodes_df)), function(i) {
      if (node_status$is_compound_node[i]) {
        # Color by compound type using status
        if (node_status$status[i] == "Input") {
          return("#7FD399")  # Green for inputs
        } else if (node_status$status[i] == "Output") {
          return("#F77370")  # Red for outputs
        } else {
          return("#B3B3B3")  # Gray for intermediates (Compound status)
        }
      } else {
        return("#BBDEFB")  # Light blue for genes
      }
    })
  }

  # Border color and style
  nodes_df$color.border <- "#333333"
  nodes_df$borderWidth <- 2

  # Dashed borders for optional nodes (not required)
  # visNetwork uses 'shapeProperties' with 'borderDashes'
  nodes_df$shapeProperties <- lapply(seq_len(nrow(nodes_df)), function(i) {
    if (!node_status$is_compound_node[i] && !node_status$required[i]) {
      list(borderDashes = c(5, 5))  # Dashed pattern: 5px dash, 5px gap
    } else {
      list(borderDashes = FALSE)  # Solid border
    }
  })

  # Build edges data frame
  edge_list <- igraph::as_edgelist(g)
  edges_df <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2],
    arrows = "to",
    color = "#999999",
    stringsAsFactors = FALSE
  )

  # Create visNetwork - full height
  vis <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "100vh")

  # Configure options
  vis <- vis %>%
    visNetwork::visNodes(
      size = 25,
      font = list(size = 14),
      borderWidthSelected = 6  # Thicker border when selected
    ) %>%
    visNetwork::visEdges(
      smooth = FALSE,  # Straight arrows
      width = 2
    ) %>%
    visNetwork::visOptions(
      manipulation = TRUE,  # Enable drag-and-drop
      highlightNearest = list(enabled = TRUE, hover = TRUE),
      nodesIdSelection = TRUE,
      autoResize = TRUE  # Adjust when window resizes
    ) %>%
    visNetwork::visInteraction(
      dragNodes = TRUE,
      dragView = TRUE,
      zoomView = TRUE,
      navigationButtons = TRUE
    ) %>%
    visNetwork::visPhysics(
      enabled = FALSE,  # Disable physics so nodes stay where placed
      stabilization = FALSE  # Don't run stabilization
    )

  # Add title
  title <- if (has_genome) {
    paste0(potato@name, " (Genome: ", attr(node_status, "genome_name"), ")")
  } else {
    potato@name
  }

  # Set layout - only set randomSeed for non-curated layouts
  if (!using_curated_coords) {
    vis <- vis %>% visNetwork::visLayout(randomSeed = 123)
  }

  # Add instructions in the title/subtitle
  title <- paste0(
    "Drag nodes to arrange, then open browser console (F12) and run: exportCoordinates()",
    if (!is.null(attr(node_status, "genome_name"))) {
      paste0(" | Genome: ", attr(node_status, "genome_name"))
    } else ""
  )

  # Store network instance and add export function
  vis <- htmlwidgets::onRender(
    vis,
    "function(el, x) {
      // Find where the network is actually stored
      console.log('el:', el);
      console.log('el.chart:', el.chart);
      console.log('this:', this);

      // Try multiple ways to access the network
      var network = el.chart || this.network || this;
      window.myVisNetwork = network;

      console.log('Stored network:', window.myVisNetwork);

      // Create global export function
      window.exportCoordinates = function() {
        console.log('exportCoordinates called');
        console.log('window.myVisNetwork:', window.myVisNetwork);

        // Try to find the network from the DOM
        var networkDiv = document.querySelector('.vis-network');
        console.log('networkDiv:', networkDiv);

        if (networkDiv && networkDiv.network) {
          var positions = networkDiv.network.getPositions();
        } else if (window.myVisNetwork && window.myVisNetwork.getPositions) {
          var positions = window.myVisNetwork.getPositions();
        } else {
          alert('Cannot find network. Check console for debug info.');
          console.error('Network not found. Available:', {
            myVisNetwork: window.myVisNetwork,
            networkDiv: networkDiv,
            hasDivNetwork: networkDiv ? networkDiv.network : null
          });
          return;
        }

        try {
          var coords = [];
          var nodeIds = Object.keys(positions);

          // Get canvas coordinates (what visNetwork is actually using internally)
          // These are absolute positions in the canvas coordinate system
          for (var i = 0; i < nodeIds.length; i++) {
            var nodeId = nodeIds[i];
            coords.push({
              id: nodeId,
              x: Math.round(positions[nodeId].x * 100) / 100,
              y: Math.round(positions[nodeId].y * 100) / 100
            });
          }

          // Download as JSON
          var dataStr = 'data:text/json;charset=utf-8,' + encodeURIComponent(JSON.stringify(coords, null, 2));
          var downloadAnchor = document.createElement('a');
          downloadAnchor.setAttribute('href', dataStr);
          downloadAnchor.setAttribute('download', 'node_coordinates.json');
          document.body.appendChild(downloadAnchor);
          downloadAnchor.click();
          downloadAnchor.remove();

          console.log('Exported ' + coords.length + ' node coordinates');
          console.log('Coordinate range - X:', Math.min(...coords.map(c => c.x)), 'to', Math.max(...coords.map(c => c.x)));
          console.log('Coordinate range - Y:', Math.min(...coords.map(c => c.y)), 'to', Math.max(...coords.map(c => c.y)));
          alert('Coordinates exported! (' + coords.length + ' nodes)');
          return coords;
        } catch (e) {
          console.error('Export error:', e);
          alert('Error: ' + e.message);
        }
      };

      console.log('Network ready. Run exportCoordinates() in console to export coordinates.');

      // Add visible export button
      setTimeout(function() {
        if (document.getElementById('export-coords-btn')) return;

        var btn = document.createElement('button');
        btn.id = 'export-coords-btn';
        btn.innerHTML = 'Export Coordinates';
        btn.style.cssText = 'position: fixed; top: 10px; right: 10px; z-index: 9999; padding: 10px 20px; background-color: #4CAF50; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 14px; box-shadow: 0 2px 5px rgba(0,0,0,0.2); font-family: Arial, sans-serif;';
        btn.onclick = function() {
          window.exportCoordinates();
        };
        document.body.appendChild(btn);
        console.log('Export button added');
      }, 1000);
    }"
  )

  vis
}


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

  # Read potato JSON
  potato_json <- jsonlite::read_json(potato_path, simplifyVector = FALSE)

  # Read coordinates
  coords <- jsonlite::read_json(coords_path, simplifyVector = TRUE)

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

  # Show which nodes will be updated
  enzyme_nodes_in_potato <- sapply(potato_json$nodes, function(n) n$id)
  enzyme_nodes_in_coords <- enzyme_nodes$id
  missing_in_coords <- setdiff(enzyme_nodes_in_potato, enzyme_nodes_in_coords)
  extra_in_coords <- setdiff(enzyme_nodes_in_coords, enzyme_nodes_in_potato)

  if (length(missing_in_coords) > 0) {
    cli::cli_alert_warning("Nodes in potato but not in coordinates: {paste(missing_in_coords, collapse=', ')}")
  }
  if (length(extra_in_coords) > 0) {
    cli::cli_alert_warning("Nodes in coordinates but not in potato: {paste(extra_in_coords, collapse=', ')}")
  }

  # Create lookup table for enzyme nodes
  enzyme_coords_lookup <- stats::setNames(
    split(enzyme_nodes[, c("x", "y")], seq_len(nrow(enzyme_nodes))),
    enzyme_nodes$id
  )

  # Determine field names based on with_compounds
  x_field <- if (with_compounds) "x_compounds" else "x"
  y_field <- if (with_compounds) "y_compounds" else "y"

  # Update enzyme nodes with coordinates
  for (i in seq_along(potato_json$nodes)) {
    node_id <- potato_json$nodes[[i]]$id

    if (node_id %in% names(enzyme_coords_lookup)) {
      potato_json$nodes[[i]][[x_field]] <- enzyme_coords_lookup[[node_id]]$x
      potato_json$nodes[[i]][[y_field]] <- enzyme_coords_lookup[[node_id]]$y
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
