#' Score pathway completeness for all genomes
#'
#' Calculates pathway completeness and confidence scores by comparing
#' annotation hits to potato pathway definitions.
#'
#' @param sack PotatoSack object with annotation results
#' @param min_fraction_override Optional. Override min_fraction threshold for all potatoes
#'
#' @returns Modified PotatoSack with scores tibble added
#'
#' @details
#' For each genome × potato combination:
#' - Maps annotation hits to potato nodes
#' - Calculates fraction of required nodes found
#' - Applies potato-specific thresholds
#' - Marks pathway as present/absent
#'
#' The scores tibble contains:
#' - genome: genome identifier
#' - potato_id: pathway identifier
#' - potato_name: pathway name
#' - nodes_found: count of detected nodes
#' - nodes_required: count of required nodes
#' - nodes_total: total nodes (required + optional)
#' - fraction: nodes_found / nodes_required
#' - present: logical, fraction >= min_fraction
#' - min_fraction: threshold used
#' - detected_genes: list column with gene IDs found
#'
#' @export
score_pathways <- function(sack, min_fraction_override = NULL) {

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("Input must be a PotatoSack object", call. = FALSE)
  }

  if (is.null(sack@results)) {
    stop("No annotation results found. Run annotate_sack() first.", call. = FALSE)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("\nScoring pathways...")
  } else {
    cli::cli_h2("Scoring pathways")
  }

  # Unnest and combine all tool results
  anno_data <- unnest_all_annotations(sack@results)

  if (nrow(anno_data) == 0) {
    warning("No annotation hits found. All pathways will be marked absent.")
  }

  # Initialize results list
  scores_list <- list()

  # For each genome
  genomes <- unique(sack@results$genome)

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("  Processing ", length(genomes), " genome(s) × ", length(sack@potatoes), " potato(es)")
  } else {
    cli::cli_alert_info("Processing {length(genomes)} genome{?s} × {length(sack@potatoes)} potato{?es}")
  }

  for (genome_id in genomes) {
    # Get hits for this genome
    genome_hits <- anno_data[anno_data$genome == genome_id, ]

    for (potato in sack@potatoes) {
      score <- score_single_pathway(
        genome_id = genome_id,
        potato = potato,
        hits = genome_hits,
        min_fraction_override = min_fraction_override
      )

      scores_list[[length(scores_list) + 1]] <- score
    }
  }

  # Combine into tibble
  scores_df <- do.call(rbind, lapply(scores_list, function(x) {
    data.frame(
      genome = x$genome,
      potato = x$potato_id,  # Use 'potato' for consistency
      potato_name = x$potato_name,
      completeness = x$fraction,  # Use 'completeness' for consistency
      nodes_found = x$nodes_found,
      nodes_required = x$nodes_required,
      nodes_total = x$nodes_total,
      present = x$present,
      present_via_marker = x$present_via_marker,
      marker_genes_found = x$marker_genes_found,
      marker_genes_total = x$marker_genes_total,
      min_fraction = x$min_fraction,
      stringsAsFactors = FALSE
    )
  }))

  # Add detected genes as list column (gene IDs, not DAG node IDs)
  scores_df$detected_genes <- lapply(scores_list, function(x) x$detected_genes)

  # Add detected marker genes as list column
  scores_df$detected_marker_genes <- lapply(scores_list, function(x) x$detected_marker_genes)

  # Add detected DAG nodes as list column
  scores_df$detected_dag_nodes <- lapply(scores_list, function(x) x$detected_dag_nodes)

  # Summary
  n_present <- sum(scores_df$present)
  n_total <- nrow(scores_df)

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("  Pathways detected: ", n_present, " / ", n_total,
            " (", round(100 * n_present / n_total, 1), "%)")
  } else {
    cli::cli_alert_success("Pathways detected: {n_present} / {n_total} ({round(100 * n_present / n_total, 1)}%)")
  }

  # Add to sack
  sack@scores <- tibble::as_tibble(scores_df)
  sack@completed_stages <- c(sack@completed_stages, "scoring")

  # Add provenance
  sack@provenance$scoring <- create_provenance(
    stage_name = "scoring",
    command = "score_pathways()",
    params = list(
      n_genomes = length(genomes),
      n_potatoes = length(sack@potatoes),
      min_fraction_override = min_fraction_override
    )
  )

  sack
}


#' Score a single genome-potato combination
#'
#' @param genome_id Character. Genome identifier
#' @param potato Potato S7 object
#' @param hits Data frame of annotation hits for this genome (all tools)
#' @param min_fraction_override Optional min_fraction override
#'
#' @returns List with scoring results
#' @keywords internal
score_single_pathway <- function(genome_id, potato, hits, min_fraction_override = NULL) {

  # Check which genes were detected and map to DAG nodes
  detected_gene_ids <- character(0)
  detected_dag_nodes <- character(0)
  detected_marker_genes <- character(0)

  for (node in potato@nodes) {
    # Get all terms for this node (gene)
    all_terms <- get_all_detection_terms_for_node(node)
    if (length(all_terms) == 0) next

    # Check if any term was hit
    gene_found <- FALSE
    for (tool_type in c("ko", "blast_term", "pfam", "hmm")) {
      term_col <- switch(tool_type,
        ko = "ko",
        blast_term = "subject",
        pfam = "pfam",
        hmm = "target",
        NULL
      )

      if (!is.null(term_col) && term_col %in% names(hits)) {
        hit_terms <- hits[[term_col]]
        if (any(all_terms %in% hit_terms)) {
          gene_found <- TRUE
          break
        }
      }
    }

    if (gene_found) {
      # Add gene ID (for reporting)
      detected_gene_ids <- c(detected_gene_ids, node$id)

      # Track marker genes separately
      is_marker <- if (!is.null(node$is_marker)) node$is_marker else FALSE
      if (is_marker) {
        detected_marker_genes <- c(detected_marker_genes, node$id)
      }

      # Add all DAG nodes for this gene (for step scoring)
      if (!is.null(node$nodes)) {
        dag_node_ids <- if (is.list(node$nodes)) unlist(node$nodes) else node$nodes
        detected_dag_nodes <- c(detected_dag_nodes, dag_node_ids)
      } else {
        # Legacy: no nodes field, use gene ID_step
        step_val <- if (is.null(node$step)) 0 else node$step
        detected_dag_nodes <- c(detected_dag_nodes, paste0(node$id, "_", step_val))
      }
    }
  }

  detected_gene_ids <- unique(detected_gene_ids)
  detected_dag_nodes <- unique(detected_dag_nodes)
  detected_marker_genes <- unique(detected_marker_genes)

  # Determine scoring method
  scoring_method <- if (!is.null(potato@scoring) && !is.null(potato@scoring$method)) {
    potato@scoring$method
  } else {
    "simple"  # Default
  }

  # Calculate scores based on method
  if (scoring_method == "step_completion") {
    # Step-based scoring: need at least one enzyme per step (uses DAG nodes)
    score_result <- score_by_steps(potato, detected_dag_nodes)
  } else {
    # Simple scoring: fraction of required genes
    score_result <- score_simple(potato, detected_gene_ids)
  }

  # Get threshold
  min_fraction <- if (!is.null(min_fraction_override)) {
    min_fraction_override
  } else if (!is.null(potato@scoring) && !is.null(potato@scoring$min_fraction)) {
    potato@scoring$min_fraction
  } else {
    0.75  # Default threshold
  }

  present <- score_result$fraction >= min_fraction

  # Check marker mode
  marker_mode <- if (!is.null(potato@scoring) && !is.null(potato@scoring$marker_mode)) {
    potato@scoring$marker_mode
  } else {
    "any"  # Default: any marker gene is sufficient
  }

  # Determine if present via marker genes
  marker_genes_total <- sum(sapply(potato@nodes, function(n) {
    if (!is.null(n$is_marker)) n$is_marker else FALSE
  }))

  # If no marker genes defined, use NA instead of 0
  if (marker_genes_total == 0) {
    present_via_marker <- NA
    marker_genes_found <- NA_integer_
    marker_genes_total <- NA_integer_
    detected_marker_genes <- character(0)  # Empty list for consistency
  } else {
    present_via_marker <- if (marker_mode == "all") {
      length(detected_marker_genes) == marker_genes_total
    } else {
      length(detected_marker_genes) > 0
    }
    marker_genes_found <- length(detected_marker_genes)
  }

  list(
    genome = genome_id,
    potato_id = potato@id,
    potato_name = potato@name,
    nodes_found = score_result$nodes_found,
    nodes_required = score_result$nodes_required,
    nodes_total = length(potato@nodes),
    fraction = score_result$fraction,
    present = present,
    present_via_marker = present_via_marker,
    marker_genes_found = marker_genes_found,
    marker_genes_total = marker_genes_total,
    min_fraction = min_fraction,
    detected_genes = detected_gene_ids,
    detected_marker_genes = detected_marker_genes,
    detected_dag_nodes = detected_dag_nodes
  )
}


#' Score pathway by step completion (handles OR branches)
#'
#' @param potato Potato S7 object
#' @param detected_dag_nodes Character vector of detected DAG node IDs (gene_step format)
#' @returns List with nodes_found, nodes_required, fraction
#' @keywords internal
score_by_steps <- function(potato, detected_dag_nodes) {

  # Extract step numbers from DAG node IDs (e.g., "aceA_3" -> "3")
  extract_step <- function(dag_node_id) {
    parts <- strsplit(dag_node_id, "_")[[1]]
    if (length(parts) >= 2) {
      return(parts[length(parts)])  # Last part is step
    }
    return(NA)
  }

  detected_steps <- unique(sapply(detected_dag_nodes, extract_step))
  detected_steps <- detected_steps[!is.na(detected_steps)]

  # Get all unique steps from potato
  all_steps <- unique(unlist(lapply(potato@nodes, function(n) {
    if (!is.null(n$step)) {
      if (is.list(n$step)) unlist(n$step) else n$step
    } else {
      NULL
    }
  })))

  all_steps <- as.character(all_steps)

  steps_complete <- sum(all_steps %in% detected_steps)
  total_steps <- length(all_steps)

  fraction <- if (total_steps > 0) steps_complete / total_steps else 0

  list(
    nodes_found = steps_complete,
    nodes_required = total_steps,
    fraction = fraction
  )
}


#' Score pathway by simple gene counting
#'
#' @param potato Potato S7 object
#' @param detected_gene_ids Character vector of detected gene IDs
#' @returns List with nodes_found, nodes_required, fraction
#' @keywords internal
score_simple <- function(potato, detected_gene_ids) {

  # Get required genes
  required_genes <- potato@nodes[sapply(potato@nodes, function(n) {
    isTRUE(n$required) || is.null(n$required)
  })]

  required_gene_ids <- sapply(required_genes, function(n) n$id)

  genes_found <- sum(required_gene_ids %in% detected_gene_ids)
  genes_required <- length(required_gene_ids)

  fraction <- if (genes_required > 0) genes_found / genes_required else 0

  list(
    nodes_found = genes_found,
    nodes_required = genes_required,
    fraction = fraction
  )
}


#' Get all detection terms for a node
#'
#' Extracts all KO/BLAST/PFAM/HMM terms from a node, handling both
#' new database schema and legacy fields.
#'
#' @param node Node list from potato
#' @returns Character vector of all terms
#' @keywords internal
get_all_detection_terms_for_node <- function(node) {
  terms <- character(0)

  # New schema: databases field
  if (!is.null(node$databases)) {
    for (db_name in names(node$databases)) {
      db_terms <- node$databases[[db_name]]
      terms <- c(terms, if (is.list(db_terms)) unlist(db_terms) else db_terms)
    }
  }

  # Legacy fields
  for (field in c("ko", "blast_terms", "pfam", "hmm")) {
    if (!is.null(node[[field]])) {
      field_terms <- node[[field]]
      terms <- c(terms, if (is.list(field_terms)) unlist(field_terms) else field_terms)
    }
  }

  unique(terms)
}
