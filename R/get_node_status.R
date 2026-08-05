#' Get node detection status for plotting (internal)
#' @noRd
get_node_status <- function(potato, sack, genome_name) {

  # Find genome index
  genome_idx <- which(sack@results$genome == genome_name)
  if (length(genome_idx) == 0) {
    cli::cli_warn("Genome {.val {genome_name}} not found in results")

    # Build graph to get node names
    g <- if (!is.null(potato@edges) && is.list(potato@edges) &&
              length(names(potato@edges)) > 0 && !is.null(potato@edges[[1]]$type)) {
      build_potato_graph(potato)
    } else {
      build_potato_graph(potato)
    }

    return(tibble::tibble(
      name = igraph::V(g)$name,
      detected = NA,
      status = "Unknown"
    ))
  }

  # Collect hits for this genome
  genome_hits <- list()

  if ("kofam" %in% names(sack@results)) {
    kofam_data <- sack@results$kofam[[genome_idx]]
    if (!is.null(kofam_data) && nrow(kofam_data) > 0) {
      # Filter by kofam thresholds
      genome_hits$kofam <- kofam_data[kofam_data$score >= kofam_data$threshold, ]
    }
  }

  if ("blast" %in% names(sack@results)) {
    blast_data <- sack@results$blast[[genome_idx]]
    if (!is.null(blast_data) && nrow(blast_data) > 0) {
      # Use default thresholds
      genome_hits$blast <- blast_data[blast_data$evalue <= 1e-10 & blast_data$bitscore >= 50, ]
    }
  }

  if ("hmm" %in% names(sack@results)) {
    hmm_data <- sack@results$hmm[[genome_idx]]
    if (!is.null(hmm_data) && nrow(hmm_data) > 0) {
      # Filter by TC or e-value
      genome_hits$hmm <- hmm_data[
        (!is.na(hmm_data$tc_threshold) & hmm_data$score >= hmm_data$tc_threshold) |
        (is.na(hmm_data$tc_threshold) & hmm_data$evalue <= 1e-10), ]
    }
  }

  # Check if multi-pathway network
  is_network <- !is.null(potato@edges) &&
                is.list(potato@edges) &&
                length(names(potato@edges)) > 0 &&
                !is.null(potato@edges[[1]]$type)

  # Get all node names from graph
  g <- build_potato_graph(potato)
  node_names <- igraph::V(g)$name

  # Check each node
  status_data <- purrr::map_dfr(node_names, function(node_name) {
    # For multi-pathway networks, node names are gene IDs directly
    # For single-pathway, extract gene ID from node_name (remove _step suffix)
    node_id <- if (is_network) node_name else sub("_\\d+$", "", node_name)

    # Find corresponding node in potato
    node <- purrr::keep(potato@nodes, ~ .x$id == node_id)

    if (length(node) == 0) {
      return(tibble::tibble(
        name = node_name,
        detected = NA,
        status = "Unknown",
        is_complex = FALSE,
        fraction_detected = NA
      ))
    }

    node <- node[[1]]

    # Check if node is detected in any tool
    detected <- FALSE

    # Check databases field
    if (!is.null(node$databases)) {
      for (db_type in names(node$databases)) {
        terms <- node$databases[[db_type]]

        if (db_type == "kofam" && !is.null(genome_hits$kofam)) {
          if (any(genome_hits$kofam$ko %in% terms)) {
            detected <- TRUE
            break
          }
        } else if (db_type == "blast" && !is.null(genome_hits$blast)) {
          if (any(genome_hits$blast$gene_id %in% terms)) {
            detected <- TRUE
            break
          }
        } else if (db_type == "hmm" && !is.null(genome_hits$hmm)) {
          if (any(genome_hits$hmm$profile %in% terms)) {
            detected <- TRUE
            break
          }
        }
      }
    }

    # Legacy: check ko field (for backward compatibility)
    if (!detected && !is.null(node$ko) && !is.null(genome_hits$kofam)) {
      if (any(genome_hits$kofam$ko %in% node$ko)) {
        detected <- TRUE
      }
    }

    # Determine status
    status <- if (detected) "Detected" else "Not detected"

    tibble::tibble(
      name = node_name,
      detected = detected,
      status = status,
      is_complex = FALSE,
      fraction_detected = if (detected) 1 else 0
    )
  })

  # Add genome name as attribute
  attr(status_data, "genome_name") <- genome_name

  status_data
}
