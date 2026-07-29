#' Score pathway presence/absence across all genomes
#'
#' Applies quality thresholds to annotation hits and calculates pathway-level
#' completion scores. Handles OR branches (alternative genes) and required vs
#' optional genes.
#'
#' @param sack PotatoSack object with annotation results
#' @param kofam_threshold Score threshold for kofam hits (NULL = use kofam threshold)
#' @param blast_evalue E-value threshold for BLAST hits (default: 1e-10)
#' @param blast_bitscore Bitscore threshold for BLAST hits (default: 50)
#' @param hmm_evalue E-value threshold for HMM hits (default: 1e-10)
#'
#' @returns Modified PotatoSack with scores in @scores
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
        # Filter by threshold
        hmm_filtered <- hmm_data[hmm_data$evalue <= hmm_evalue, ]
        genome_hits$hmm <- hmm_filtered
      }
    }

    # Score each potato for this genome
    for (potato in sack@potatoes) {
      potato_score <- score_single_pathway(
        potato = potato,
        genome_name = genome_name,
        genome_hits = genome_hits
      )

      all_scores[[length(all_scores) + 1]] <- potato_score
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
    if (is.list(n$step)) n$step[[1]] else n$step
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

  # Calculate completion
  steps_detected <- sum(unlist(step_completion))
  steps_total <- length(steps)
  fraction <- steps_detected / steps_total

  # Determine presence based on min_fraction threshold
  min_fraction <- potato@scoring$min_fraction
  if (is.null(min_fraction)) min_fraction <- 0.75

  present <- fraction >= min_fraction

  # Return score
  list(
    genome = genome_name,
    potato = potato@id,
    potato_name = potato@name,
    steps_detected = steps_detected,
    steps_total = steps_total,
    fraction = fraction,
    present = present
  )
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
