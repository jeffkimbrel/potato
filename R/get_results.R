#' Get gene-level annotation results
#'
#' Returns a tibble of all gene-level annotation hits across all tools,
#' formatted as a long table compatible with GATOR output format. Includes
#' a "passed" column indicating whether each hit passes quality thresholds.
#'
#' @param sack PotatoSack object with annotation results
#'
#' @return Tibble with columns: genome, potato, node_id, gene_id, tool, passed, and tool-specific columns
#' @export
get_gene_results <- function(sack) {

  if (!S7::S7_inherits(sack, PotatoSack)) {
    cli::cli_abort("Input must be a PotatoSack object")
  }

  if (is.null(sack@results) || nrow(sack@results) == 0) {
    cli::cli_abort("No annotation results found in sack. Run annotation tools first (run_kofam, run_blast, run_hmm).")
  }

  # Get threshold settings from config
  blast_evalue_threshold <- sack@config$databases$blast$thresholds$evalue %||% 1e-10
  blast_bitscore_threshold <- sack@config$databases$blast$thresholds$bitscore %||% 50
  hmm_evalue_threshold <- sack@config$databases$hmm$thresholds$evalue %||% 1e-10

  # Combine results from all tools
  all_results <- list()

  # KOfam results (has threshold column per gene)
  if ("kofam" %in% names(sack@results) && any(!sapply(sack@results$kofam, is.null))) {
    kofam_data <- sack@results %>%
      dplyr::select(genome, kofam) %>%
      tidyr::unnest(cols = kofam) %>%
      dplyr::mutate(
        tool = "kofam",
        passed = score >= threshold
      )

    all_results$kofam <- kofam_data
  }

  # BLAST results (uses global thresholds)
  if ("blast" %in% names(sack@results) && any(!sapply(sack@results$blast, is.null))) {
    blast_data <- sack@results %>%
      dplyr::select(genome, blast) %>%
      tidyr::unnest(cols = blast) %>%
      dplyr::mutate(
        tool = "blast",
        passed = evalue <= blast_evalue_threshold & bitscore >= blast_bitscore_threshold
      )

    all_results$blast <- blast_data
  }

  # HMM results (uses tc_threshold if available, else evalue)
  if ("hmm" %in% names(sack@results) && any(!sapply(sack@results$hmm, is.null))) {
    hmm_data <- sack@results %>%
      dplyr::select(genome, hmm) %>%
      tidyr::unnest(cols = hmm) %>%
      dplyr::mutate(
        tool = "hmm",
        passed = dplyr::case_when(
          !is.na(tc_threshold) ~ score >= tc_threshold,  # Use TC if available
          TRUE ~ evalue <= hmm_evalue_threshold          # Otherwise use evalue
        )
      )

    all_results$hmm <- hmm_data
  }

  if (length(all_results) == 0) {
    cli::cli_abort("No annotation data to export. Run annotation tools first (run_kofam, run_blast, run_hmm).")
  }

  # Combine all tools - use bind_rows which handles different columns
  combined <- dplyr::bind_rows(all_results)

  # Reorder columns to put key ones first (including passed)
  key_cols <- c("genome", "potato", "node_id", "gene_id", "tool", "passed")
  other_cols <- setdiff(names(combined), key_cols)
  combined <- combined %>%
    dplyr::select(dplyr::all_of(key_cols[key_cols %in% names(combined)]), dplyr::all_of(other_cols))

  combined
}


#' Get pathway-level scoring results
#'
#' Returns a tibble of pathway completion scores across all genomes,
#' formatted for analysis and export. Includes potato_hash for version tracking.
#'
#' @param sack PotatoSack object with scoring results
#'
#' @return Tibble with columns: genome, potato, potato_hash, pathway_name, total_steps_detected,
#'   total_steps, fraction, min_fraction, present, essential_total_steps_detected,
#'   essential_steps, essential_fraction, essential_pathway_present
#' @export
get_pathway_scores <- function(sack) {

  if (!S7::S7_inherits(sack, PotatoSack)) {
    cli::cli_abort("Input must be a PotatoSack object")
  }

  if (is.null(sack@scores) || nrow(sack@scores) == 0) {
    cli::cli_abort("No scoring results found in sack. Run score_pathways() first.")
  }

  scores <- sack@scores

  # Add potato_hash from potatoes if not already present
  if (!"potato_hash" %in% names(scores)) {
    potato_hashes <- sapply(sack@potatoes, function(p) compute_potato_hash(p))
    names(potato_hashes) <- names(sack@potatoes)

    scores <- scores %>%
      dplyr::mutate(potato_hash = potato_hashes[potato])
  }

  # Reorder columns to put potato_hash early
  key_cols <- c("genome", "potato", "potato_hash", "pathway_name")
  other_cols <- setdiff(names(scores), key_cols)
  scores <- scores %>%
    dplyr::select(dplyr::all_of(key_cols[key_cols %in% names(scores)]), dplyr::all_of(other_cols))

  # Remove pathway and potato columns if they exist
  if ("pathway" %in% names(scores)) {
    scores <- scores %>% dplyr::select(-pathway)
  }
  if ("potato" %in% names(scores)) {
    scores <- scores %>% dplyr::select(-potato)
  }

  scores
}
