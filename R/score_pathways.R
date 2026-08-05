#' Score pathway presence/absence across all genomes
#'
#' Applies quality thresholds to annotation hits and calculates pathway-level
#' completion scores. For multi-pathway networks, scores each pathway independently.
#'
#' Handles OR branches (alternative genes), required vs optional genes, and
#' multi-pathway networks where genes are shared across pathways.
#'
#' Threshold priority:
#' - Kofam: Uses per-gene threshold from KEGG (can override with kofam_threshold)
#' - HMM: Uses per-profile TC (trusted cutoff) if available, otherwise hmm_evalue
#' - BLAST: Uses global blast_evalue and blast_bitscore
#'
#' @param sack PotatoSack object with annotation results
#' @param kofam_threshold Score threshold for kofam hits (NULL = use per-gene threshold)
#' @param blast_evalue E-value threshold for BLAST hits (default: 1e-10)
#' @param blast_bitscore Bitscore threshold for BLAST hits (default: 50)
#' @param hmm_evalue E-value threshold for HMM hits without TC (default: 1e-10)
#'
#' @returns Modified PotatoSack with scores in @scores. For multi-pathway networks,
#'   scores tibble includes 'pathway' and 'pathway_name' columns with one row per
#'   pathway per genome. Scoring includes both all-steps metrics (total_steps_detected,
#'   total_steps, fraction, present) and required-only metrics (essential_total_steps_detected,
#'   essential_steps, essential_fraction, essential_pathway_present).
#' @export

score_pathways <- function(sack,
                           kofam_threshold = NULL,
                           blast_evalue = 1e-10,
                           blast_bitscore = 50,
                           hmm_evalue = 1e-10) {

  # Check that we have annotation results
  if (is.null(sack@results)) {
    cli::cli_abort("No annotation results found. Run annotation tools first.")
  }

  cli::cli_alert_info("Scoring pathways across {nrow(sack@results)} genome{?s}...")

  # Get genome names
  genome_names <- sack@results$genome

  # Initialize scores list
  all_scores <- list()

  # Process each genome
  for (i in seq_along(genome_names)) {
    genome_name <- genome_names[i]

    # Collect all hits for this genome across all tools
    genome_hits <- list()

    # Kofam hits
    if ("kofam" %in% names(sack@results)) {
      kofam_data <- sack@results$kofam[[i]]
      if (!is.null(kofam_data) && nrow(kofam_data) > 0) {
        # Filter by threshold
        if (is.null(kofam_threshold)) {
          # Use kofam's own threshold
          kofam_filtered <- kofam_data[kofam_data$score >= kofam_data$threshold, ]
        } else {
          kofam_filtered <- kofam_data[kofam_data$score >= kofam_threshold, ]
        }
        genome_hits$kofam <- kofam_filtered
      }
    }

    # BLAST hits
    if ("blast" %in% names(sack@results)) {
      blast_data <- sack@results$blast[[i]]
      if (!is.null(blast_data) && nrow(blast_data) > 0) {
        # Filter by thresholds
        blast_filtered <- blast_data[
          blast_data$evalue <= blast_evalue &
          blast_data$bitscore >= blast_bitscore, ]
        genome_hits$blast <- blast_filtered
      }
    }

    # HMM hits
    if ("hmm" %in% names(sack@results)) {
      hmm_data <- sack@results$hmm[[i]]
      if (!is.null(hmm_data) && nrow(hmm_data) > 0) {
        # Filter by threshold: use TC if available, otherwise use e-value
        hmm_filtered <- hmm_data[
          (!is.na(hmm_data$tc_threshold) & hmm_data$score >= hmm_data$tc_threshold) |
          (is.na(hmm_data$tc_threshold) & hmm_data$evalue <= hmm_evalue), ]
        genome_hits$hmm <- hmm_filtered
      }
    }

    # Score each potato for this genome
    for (potato in sack@potatoes) {
      # Check if this is a multi-pathway network
      is_network <- !is.null(potato@edges) &&
                    is.list(potato@edges) &&
                    length(names(potato@edges)) > 0 &&
                    !is.null(potato@edges[[1]]$type)  # Pathway has type field

      if (is_network) {
        # Score each pathway independently
        for (pathway_id in names(potato@edges)) {
          pathway <- potato@edges[[pathway_id]]

          pathway_score <- score_single_pathway_network(
            potato_id = potato@id,
            potato_name = potato@name,
            pathway_id = pathway_id,
            pathway = pathway,
            global_nodes = potato@nodes,
            genome_name = genome_name,
            genome_hits = genome_hits
          )

          all_scores[[length(all_scores) + 1]] <- pathway_score
        }
      } else {
        # Single-pathway potato (legacy schema)
        potato_score <- score_single_pathway(
          potato = potato,
          genome_name = genome_name,
          genome_hits = genome_hits
        )

        all_scores[[length(all_scores) + 1]] <- potato_score
      }
    }
  }

  # Combine into tibble
  sack@scores <- dplyr::bind_rows(all_scores)

  cli::cli_alert_success("Scored {length(sack@potatoes)} pathway{?s} across {length(genome_names)} genome{?s}")

  sack
}


#' Score a single pathway for a single genome (internal)
#' @noRd
score_single_pathway <- function(potato, genome_name, genome_hits) {

  # Get all nodes in pathway
  nodes <- potato@nodes

  # Group nodes by step (for handling OR branches)
  steps <- unique(sapply(nodes, function(n) {
    step_val <- n$step
    if (is.list(step_val)) {
      # If step is a list, get first element
      step_val <- step_val[[1]]
    }
    # Ensure it's numeric
    as.numeric(step_val)
  }))
  steps <- sort(steps)

  # Track which steps are complete
  step_completion <- list()

  for (step_num in steps) {
    # Get all nodes at this step (OR alternatives)
    step_nodes <- Filter(function(n) {
      s <- if (is.list(n$step)) n$step[[1]] else n$step
      s == step_num
    }, nodes)

    # Check if ANY node at this step is detected
    step_detected <- FALSE

    for (node in step_nodes) {
      node_detected <- is_node_detected(node, potato@id, genome_hits)
      if (node_detected) {
        step_detected <- TRUE
        break  # OR branch satisfied
      }
    }

    step_completion[[as.character(step_num)]] <- step_detected
  }

  # Calculate completion (all steps)
  total_steps_detected <- sum(unlist(step_completion))
  total_steps <- length(steps)
  fraction <- total_steps_detected / total_steps

  # Determine presence based on min_fraction threshold
  min_fraction <- potato@scoring$min_fraction
  if (is.null(min_fraction)) min_fraction <- 0.75

  present <- fraction >= min_fraction

  # Calculate completion for required steps only
  required_steps <- Filter(function(step_num) {
    step_nodes <- Filter(function(n) {
      s <- if (is.list(n$step)) n$step[[1]] else n$step
      s == step_num
    }, nodes)
    # Step is required if ANY node at that step is required
    any(sapply(step_nodes, function(n) n$required %||% FALSE))
  }, steps)

  if (length(required_steps) > 0) {
    essential_total_steps_detected <- sum(unlist(step_completion[as.character(required_steps)]))
    essential_steps <- length(required_steps)
    essential_fraction <- essential_total_steps_detected / essential_steps
    essential_pathway_present <- essential_fraction >= min_fraction
  } else {
    # No required steps defined
    essential_total_steps_detected <- NA_integer_
    essential_steps <- NA_integer_
    essential_fraction <- NA_real_
    essential_pathway_present <- NA
  }

  # Return score
  list(
    genome = genome_name,
    potato = potato@id,
    potato_name = potato@name,
    total_steps_detected = total_steps_detected,
    total_steps = total_steps,
    fraction = fraction,
    min_fraction = min_fraction,
    present = present,
    essential_total_steps_detected = essential_total_steps_detected,
    essential_steps = essential_steps,
    essential_fraction = essential_fraction,
    essential_pathway_present = essential_pathway_present
  )
}


#' Score a single pathway in a multi-pathway network (internal)
#' @noRd
score_single_pathway_network <- function(potato_id, potato_name, pathway_id,
                                         pathway, global_nodes, genome_name,
                                         genome_hits) {

  # Get pathway-specific nodes
  pathway_nodes <- pathway$nodes

  if (is.null(pathway_nodes) || length(pathway_nodes) == 0) {
    # Empty pathway
    return(list(
      genome = genome_name,
      potato = potato_id,
      potato_name = potato_name,
      pathway = pathway_id,
      pathway_name = pathway$name,
      total_steps_detected = 0,
      total_steps = 0,
      fraction = 0,
      present = FALSE
    ))
  }

  # Build merged nodes: global detection methods + pathway-specific attributes
  merged_nodes <- list()
  for (node_id in names(pathway_nodes)) {
    # Find global node
    global_node <- Find(function(n) n$id == node_id, global_nodes)

    if (is.null(global_node)) {
      cli::cli_warn("Pathway '{pathway_id}' references node '{node_id}' not found in global nodes")
      next
    }

    # Merge: global databases + pathway-specific step/required/marker
    merged_node <- global_node
    merged_node$step <- pathway_nodes[[node_id]]$step
    merged_node$required <- pathway_nodes[[node_id]]$required
    merged_node$marker <- pathway_nodes[[node_id]]$marker

    merged_nodes[[length(merged_nodes) + 1]] <- merged_node
  }

  # Group nodes by step (for handling OR branches)
  steps <- unique(sapply(merged_nodes, function(n) {
    step_val <- n$step
    if (is.list(step_val)) {
      # If step is a list, get first element
      step_val <- step_val[[1]]
    }
    # Ensure it's numeric
    as.numeric(step_val)
  }))
  steps <- sort(steps)

  # Track which steps are complete
  step_completion <- list()

  for (step_num in steps) {
    # Get all nodes at this step (OR alternatives)
    step_nodes <- Filter(function(n) {
      s <- if (is.list(n$step)) n$step[[1]] else n$step
      s == step_num
    }, merged_nodes)

    # Check if ANY node at this step is detected
    step_detected <- FALSE

    for (node in step_nodes) {
      node_detected <- is_node_detected_network(node, potato_id, genome_hits)
      if (node_detected) {
        step_detected <- TRUE
        break  # OR branch satisfied
      }
    }

    step_completion[[as.character(step_num)]] <- step_detected
  }

  # Calculate completion (all steps)
  total_steps_detected <- sum(unlist(step_completion))
  total_steps <- length(steps)
  fraction <- total_steps_detected / total_steps

  # Determine presence based on min_fraction threshold
  min_fraction <- pathway$scoring$min_fraction
  if (is.null(min_fraction)) min_fraction <- 0.75

  present <- fraction >= min_fraction

  # Calculate completion for required steps only
  required_steps <- Filter(function(step_num) {
    step_nodes <- Filter(function(n) {
      s <- if (is.list(n$step)) n$step[[1]] else n$step
      s == step_num
    }, merged_nodes)
    # Step is required if ANY node at that step is required
    any(sapply(step_nodes, function(n) n$required %||% FALSE))
  }, steps)

  if (length(required_steps) > 0) {
    essential_total_steps_detected <- sum(unlist(step_completion[as.character(required_steps)]))
    essential_steps <- length(required_steps)
    essential_fraction <- essential_total_steps_detected / essential_steps
    essential_pathway_present <- essential_fraction >= min_fraction
  } else {
    # No required steps defined
    essential_total_steps_detected <- NA_integer_
    essential_steps <- NA_integer_
    essential_fraction <- NA_real_
    essential_pathway_present <- NA
  }

  # Return score
  list(
    genome = genome_name,
    potato = potato_id,
    potato_name = potato_name,
    pathway = pathway_id,
    pathway_name = pathway$name %||% pathway_id,
    total_steps_detected = total_steps_detected,
    total_steps = total_steps,
    fraction = fraction,
    min_fraction = min_fraction,
    present = present,
    essential_total_steps_detected = essential_total_steps_detected,
    essential_steps = essential_steps,
    essential_fraction = essential_fraction,
    essential_pathway_present = essential_pathway_present
  )
}


#' Check if a node is detected in genome hits for network pathways (internal)
#' @noRd
is_node_detected_network <- function(node, potato_id, genome_hits) {

  # Check each database type
  databases <- node$databases

  if (is.null(databases)) return(FALSE)

  # For network pathways, hits are stored with potato_id but we need to match
  # by node_id (gene ID) since genes are shared across pathways

  # Check kofam
  if (!is.null(databases$kofam) && !is.null(genome_hits$kofam)) {
    kofam_hits <- genome_hits$kofam
    # Match by potato_id and node_id
    node_hits <- kofam_hits[
      kofam_hits$potato == potato_id &
      kofam_hits$node_id == node$id, ]
    if (nrow(node_hits) > 0) return(TRUE)
  }

  # Check blast
  if (!is.null(databases$blast) && !is.null(genome_hits$blast)) {
    blast_hits <- genome_hits$blast
    node_hits <- blast_hits[
      blast_hits$potato == potato_id &
      blast_hits$node_id == node$id, ]
    if (nrow(node_hits) > 0) return(TRUE)
  }

  # Check hmm
  if (!is.null(databases$hmm) && !is.null(genome_hits$hmm)) {
    hmm_hits <- genome_hits$hmm
    node_hits <- hmm_hits[
      hmm_hits$potato == potato_id &
      hmm_hits$node_id == node$id, ]
    if (nrow(node_hits) > 0) return(TRUE)
  }

  return(FALSE)
}


#' Check if a node is detected in genome hits (internal)
#' @noRd
is_node_detected <- function(node, potato_id, genome_hits) {

  # Check each database type
  databases <- node$databases

  if (is.null(databases)) return(FALSE)

  # Check kofam
  if (!is.null(databases$kofam) && !is.null(genome_hits$kofam)) {
    kofam_hits <- genome_hits$kofam
    # Check if this node appears in filtered kofam hits
    node_hits <- kofam_hits[
      kofam_hits$potato == potato_id &
      kofam_hits$node_id == node$id, ]
    if (nrow(node_hits) > 0) return(TRUE)
  }

  # Check blast
  if (!is.null(databases$blast) && !is.null(genome_hits$blast)) {
    blast_hits <- genome_hits$blast
    node_hits <- blast_hits[
      blast_hits$potato == potato_id &
      blast_hits$node_id == node$id, ]
    if (nrow(node_hits) > 0) return(TRUE)
  }

  # Check hmm
  if (!is.null(databases$hmm) && !is.null(genome_hits$hmm)) {
    hmm_hits <- genome_hits$hmm
    node_hits <- hmm_hits[
      hmm_hits$potato == potato_id &
      hmm_hits$node_id == node$id, ]
    if (nrow(node_hits) > 0) return(TRUE)
  }

  return(FALSE)
}
