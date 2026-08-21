#' Identify genes missing across genomes
#'
#' Analyzes which genes are systematically missing, with distinction between
#' genes with no hits, hits below threshold, and hits passing threshold.
#'
#' @param sack PotatoSack object with annotation results and scores
#' @param potato_name Optional: filter to specific pathway
#' @param min_genomes Minimum genomes to consider (default: 3)
#'
#' @returns Tibble showing genes and detection status across thresholds
#' @export
summarize_missing_genes <- function(sack, potato_name = NULL, min_genomes = 3) {

  if (is.null(sack@results)) {
    cli::cli_abort("No annotation results found. Run annotation tools first.")
  }

  # Filter to specific potato if requested
  potatoes <- sack@potatoes
  if (!is.null(potato_name)) {
    potatoes <- potatoes[sapply(potatoes, function(p) p@name == potato_name)]
  }

  # Helper: Check if gene has any annotation hits
  check_gene_has_hits <- function(gene, potato_id, genome_hits) {
    databases <- gene$databases
    if (is.null(databases)) return(FALSE)

    # Check kofam
    if (!is.null(databases$kofam) && !is.null(genome_hits$kofam)) {
      gene_hits <- genome_hits$kofam[
        genome_hits$kofam$potato == potato_id &
        genome_hits$kofam$node_id == gene$id, ]
      if (nrow(gene_hits) > 0) return(TRUE)
    }

    # Check blast
    if (!is.null(databases$blast) && !is.null(genome_hits$blast)) {
      gene_hits <- genome_hits$blast[
        genome_hits$blast$potato == potato_id &
        genome_hits$blast$node_id == gene$id, ]
      if (nrow(gene_hits) > 0) return(TRUE)
    }

    # Check hmm
    if (!is.null(databases$hmm) && !is.null(genome_hits$hmm)) {
      gene_hits <- genome_hits$hmm[
        genome_hits$hmm$potato == potato_id &
        genome_hits$hmm$node_id == gene$id, ]
      if (nrow(gene_hits) > 0) return(TRUE)
    }

    return(FALSE)
  }

  # Helper: Check if gene has passing hits (applies thresholds)
  check_gene_passes_threshold <- function(gene, potato_id, genome_hits) {
    databases <- gene$databases
    if (is.null(databases)) return(FALSE)

    # Check kofam (score >= threshold)
    if (!is.null(databases$kofam) && !is.null(genome_hits$kofam)) {
      gene_hits <- genome_hits$kofam[
        genome_hits$kofam$potato == potato_id &
        genome_hits$kofam$node_id == gene$id, ]
      if (nrow(gene_hits) > 0) {
        if (any(gene_hits$score >= gene_hits$threshold)) return(TRUE)
      }
    }

    # Check blast (evalue <= 1e-10 AND bitscore >= 50)
    if (!is.null(databases$blast) && !is.null(genome_hits$blast)) {
      gene_hits <- genome_hits$blast[
        genome_hits$blast$potato == potato_id &
        genome_hits$blast$node_id == gene$id, ]
      if (nrow(gene_hits) > 0) {
        if (any(gene_hits$evalue <= 1e-10 & gene_hits$bitscore >= 50)) return(TRUE)
      }
    }

    # Check hmm (score >= TC or evalue <= 1e-10)
    if (!is.null(databases$hmm) && !is.null(genome_hits$hmm)) {
      gene_hits <- genome_hits$hmm[
        genome_hits$hmm$potato == potato_id &
        genome_hits$hmm$node_id == gene$id, ]
      if (nrow(gene_hits) > 0) {
        passing <- (!is.na(gene_hits$tc_threshold) & gene_hits$score >= gene_hits$tc_threshold) |
                   (is.na(gene_hits$tc_threshold) & gene_hits$evalue <= 1e-10)
        if (any(passing)) return(TRUE)
      }
    }

    return(FALSE)
  }

  # Track gene detection status
  detection_data <- list()

  for (genome_name in sack@results$genome) {
    genome_idx <- which(sack@results$genome == genome_name)

    # Collect hits (raw data with thresholds)
    genome_hits_raw <- list()
    if ("kofam" %in% names(sack@results)) genome_hits_raw$kofam <- sack@results$kofam[[genome_idx]]
    if ("blast" %in% names(sack@results)) genome_hits_raw$blast <- sack@results$blast[[genome_idx]]
    if ("hmm" %in% names(sack@results)) genome_hits_raw$hmm <- sack@results$hmm[[genome_idx]]

    # Check each potato
    for (potato in potatoes) {
      for (gene in potato@genes) {
        # Check if gene has ANY hits (passing or not)
        has_any_hits <- check_gene_has_hits(gene, potato@id, genome_hits_raw)

        # Check if gene has PASSING hits (apply thresholds)
        has_passing_hits <- check_gene_passes_threshold(gene, potato@id, genome_hits_raw)

        detection_data[[length(detection_data) + 1]] <- list(
          potato = potato@name,
          gene_id = gene$id,
          gene_name = gene$name,
          genome = genome_name,
          has_any_hits = has_any_hits,
          has_passing_hits = has_passing_hits
        )
      }
    }
  }

  # Convert to tibble and summarize
  result <- dplyr::bind_rows(detection_data) %>%
    dplyr::group_by(potato, gene_id, gene_name) %>%
    dplyr::summarize(
      total_genomes = dplyr::n(),
      times_no_hits = sum(!has_any_hits),
      times_below_threshold = sum(has_any_hits & !has_passing_hits),
      times_passing = sum(has_passing_hits),
      fraction_passing = times_passing / total_genomes,
      .groups = "drop"
    ) %>%
    dplyr::filter(total_genomes >= min_genomes) %>%
    dplyr::arrange(fraction_passing)

  result
}


#' Find pathways that are close to being called present
#'
#' Identifies "near miss" pathways that are just below the threshold,
#' which may indicate threshold tuning issues.
#'
#' @param sack PotatoSack object with scores
#' @param buffer Distance from threshold to consider "near" (default: 0.1)
#'
#' @returns Tibble of near-miss pathways
#' @export
find_near_miss_pathways <- function(sack, buffer = 0.1) {

  if (is.null(sack@scores)) {
    cli::cli_abort("No scores found. Run {.fn score_pathways} first.")
  }

  # Find pathways that are close but not present
  near_miss <- sack@scores %>%
    dplyr::filter(!present) %>%
    dplyr::mutate(
      # Use min_fraction already in scores tibble (added during scoring)
      distance_from_threshold = min_fraction - fraction,
      near_miss = distance_from_threshold <= buffer
    ) %>%
    dplyr::filter(near_miss) %>%
    dplyr::select(genome, potato_name, pathway_name, fraction, min_fraction,
                  distance_from_threshold, total_steps_detected, total_steps) %>%
    dplyr::arrange(distance_from_threshold)

  if (nrow(near_miss) == 0) {
    cli::cli_alert_info("No near-miss pathways found (within {buffer} of threshold)")
  } else {
    cli::cli_alert_success("Found {nrow(near_miss)} near-miss pathway{?s}")
  }

  near_miss
}



#' Inspect gene-level annotation hits and thresholds
#'
#' Shows all annotation hits for genes in a pathway, with scores vs thresholds.
#' Helps identify which genes are failing thresholds and by how much.
#'
#' @param sack PotatoSack object with annotation results
#' @param potato_name Optional: filter to specific pathway
#' @param show_passing Include hits that pass thresholds (default: FALSE)
#'
#' @returns Tibble with columns: potato, gene_id, gene_name, genome, tool,
#'   score, threshold, margin, passed
#' @export
inspect_gene_thresholds <- function(sack, potato_name = NULL, show_passing = FALSE) {

  if (is.null(sack@results)) {
    cli::cli_abort("No annotation results found. Run annotation tools first.")
  }

  # Filter to specific potato if requested
  potatoes <- sack@potatoes
  if (!is.null(potato_name)) {
    potatoes <- potatoes[sapply(potatoes, function(p) p@name == potato_name)]
    if (length(potatoes) == 0) {
      cli::cli_abort("Potato {.val {potato_name}} not found")
    }
  }

  # Collect all hits with threshold info
  hit_data <- list()

  for (genome_name in sack@results$genome) {
    genome_idx <- which(sack@results$genome == genome_name)

    # Check each potato
    for (potato in potatoes) {
      for (gene in potato@genes) {

        # Check kofam
        if (!is.null(gene$databases$kofam) && "kofam" %in% names(sack@results)) {
          kofam_hits <- sack@results$kofam[[genome_idx]]
          if (!is.null(kofam_hits) && nrow(kofam_hits) > 0) {
            gene_hits <- kofam_hits[
              kofam_hits$potato == potato@id &
              kofam_hits$node_id == gene$id, ]

            if (nrow(gene_hits) > 0) {
              for (i in 1:nrow(gene_hits)) {
                hit <- gene_hits[i, ]
                passed <- hit$score >= hit$threshold
                hit_data[[length(hit_data) + 1]] <- list(
                  potato = potato@name,
                  gene_id = gene$id,
                  gene_name = gene$name,
                  genome = genome_name,
                  tool = "kofam",
                  score = hit$score,
                  threshold = hit$threshold,
                  margin = hit$score - hit$threshold,
                  passed = passed
                )
              }
            }
          }
        }

        # Check BLAST
        if (!is.null(gene$databases$blast) && "blast" %in% names(sack@results)) {
          blast_hits <- sack@results$blast[[genome_idx]]
          if (!is.null(blast_hits) && nrow(blast_hits) > 0) {
            gene_hits <- blast_hits[
              blast_hits$potato == potato@id &
              blast_hits$node_id == gene$id, ]

            if (nrow(gene_hits) > 0) {
              # BLAST uses global thresholds (from scoring default: 1e-10, 50)
              # Show evalue as "score" for consistency
              for (i in 1:nrow(gene_hits)) {
                hit <- gene_hits[i, ]
                passed <- hit$evalue <= 1e-10 & hit$bitscore >= 50
                hit_data[[length(hit_data) + 1]] <- list(
                  potato = potato@name,
                  gene_id = gene$id,
                  gene_name = gene$name,
                  genome = genome_name,
                  tool = "blast",
                  score = hit$evalue,  # Show evalue as score
                  threshold = 1e-10,   # Default threshold
                  margin = 1e-10 - hit$evalue,  # Negative = passing for evalue
                  passed = passed
                )
              }
            }
          }
        }

        # Check HMM
        if (!is.null(gene$databases$hmm) && "hmm" %in% names(sack@results)) {
          hmm_hits <- sack@results$hmm[[genome_idx]]
          if (!is.null(hmm_hits) && nrow(hmm_hits) > 0) {
            gene_hits <- hmm_hits[
              hmm_hits$potato == potato@id &
              hmm_hits$node_id == gene$id, ]

            if (nrow(gene_hits) > 0) {
              for (i in 1:nrow(gene_hits)) {
                hit <- gene_hits[i, ]
                # HMM uses TC if available, else evalue
                if (!is.na(hit$tc_threshold)) {
                  passed <- hit$score >= hit$tc_threshold
                  threshold <- hit$tc_threshold
                  margin <- hit$score - hit$tc_threshold
                } else {
                  passed <- hit$evalue <= 1e-10
                  threshold <- 1e-10
                  margin <- 1e-10 - hit$evalue  # Negative = passing for evalue
                }

                hit_data[[length(hit_data) + 1]] <- list(
                  potato = potato@name,
                  gene_id = gene$id,
                  gene_name = gene$name,
                  genome = genome_name,
                  tool = "hmm",
                  score = if (!is.na(hit$tc_threshold)) hit$score else hit$evalue,
                  threshold = threshold,
                  margin = margin,
                  passed = passed
                )
              }
            }
          }
        }
      }
    }
  }

  # Convert to tibble
  if (length(hit_data) == 0) {
    cli::cli_alert_warning("No annotation hits found")
    return(tibble::tibble())
  }

  result <- dplyr::bind_rows(hit_data)

  # If not showing passing hits, only show failing hits for genes with NO passing hits
  if (!show_passing) {
    # Find genome + gene combos that have at least one passing hit
    genes_with_passing <- result %>%
      dplyr::filter(passed) %>%
      dplyr::distinct(potato, gene_id, genome)

    # Remove failing hits for genes that already have a passing hit
    result <- result %>%
      dplyr::anti_join(genes_with_passing, by = c("potato", "gene_id", "genome")) %>%
      dplyr::filter(!passed)

    # For remaining failing hits, keep only the best (closest to passing) per gene/genome
    result <- result %>%
      dplyr::group_by(potato, gene_id, genome) %>%
      dplyr::slice_min(abs(margin), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  # Sort by absolute margin (closest to passing first)
  result <- result %>% dplyr::arrange(abs(margin))

  # Summary message
  n_total <- length(hit_data)
  n_failing <- sum(!result$passed)
  n_passing <- sum(result$passed)

  if (show_passing) {
    cli::cli_alert_info("Found {n_total} hit{?s}: {n_passing} passing, {n_failing} failing")
  } else {
    cli::cli_alert_info("Found {n_failing} failing hit{?s} (use show_passing = TRUE to see all)")
  }

  result
}
