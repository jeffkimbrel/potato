#' Print provenance summary
#'
#' Displays provenance information for all annotation and scoring steps
#'
#' @param sack PotatoSack object
#' @export
print_provenance <- function(sack) {

  if (is.null(sack@provenance) || length(sack@provenance) == 0) {
    cli::cli_alert_info("No provenance data found")
    return(invisible(NULL))
  }

  cli::cli_h1("POTATO Provenance")

  # Genomes
  if (!is.null(sack@provenance$genomes)) {
    cli::cli_h2("Genomes Added")
    for (i in seq_along(sack@provenance$genomes)) {
      genome_prov <- sack@provenance$genomes[[i]]
      cli::cli_text("{i}. {genome_prov$timestamp}: Added {genome_prov$n_added} genome{?s}")
      if (length(genome_prov$genome_names) <= 5) {
        cli::cli_ul(genome_prov$genome_names)
      } else {
        cli::cli_text("   First 5: {paste(head(genome_prov$genome_names, 5), collapse = ', ')}...")
      }
    }
    cli::cli_text("")
  }

  # Kofam
  if (!is.null(sack@provenance$kofam)) {
    kofam <- sack@provenance$kofam
    cli::cli_h2("Kofam Annotation")
    cli::cli_dl(c(
      "Timestamp" = kofam$timestamp,
      "Tool version" = kofam$tool_version,
      "Conda env" = kofam$conda_env %||% "none",
      "Workers" = as.character(kofam$workers),
      "Potatoes requested" = paste(kofam$potatoes_requested, collapse = ", "),
      "Potatoes with KOs" = paste(kofam$potatoes_with_genes, collapse = ", "),
      "Genomes" = as.character(kofam$n_genomes),
      "KO terms" = as.character(kofam$n_kos)
    ))
    cli::cli_text("")
    cli::cli_text("Command template:")
    cli::cli_code(kofam$commands$exec_annotation_template)
    cli::cli_text("")
  }

  # BLAST
  if (!is.null(sack@provenance$blast)) {
    blast <- sack@provenance$blast
    cli::cli_h2("BLAST Annotation")
    cli::cli_dl(c(
      "Timestamp" = blast$timestamp,
      "Tool version" = blast$tool_version,
      "Conda env" = blast$conda_env %||% "none",
      "Workers" = as.character(blast$workers),
      "Potatoes requested" = paste(blast$potatoes_requested, collapse = ", "),
      "Potatoes with BLAST" = paste(blast$potatoes_with_genes, collapse = ", "),
      "Genomes" = as.character(blast$n_genomes),
      "Subject sequences" = as.character(blast$n_subjects)
    ))
    cli::cli_text("")
    cli::cli_text("Command templates:")
    cli::cli_code(blast$commands$makeblastdb_template)
    cli::cli_code(blast$commands$blastp_template)
    cli::cli_text("")
  }

  # HMM
  if (!is.null(sack@provenance$hmm)) {
    hmm <- sack@provenance$hmm
    cli::cli_h2("HMM Annotation")
    cli::cli_dl(c(
      "Timestamp" = hmm$timestamp,
      "Tool version" = hmm$tool_version,
      "Conda env" = hmm$conda_env %||% "none",
      "Workers" = as.character(hmm$workers),
      "Potatoes requested" = paste(hmm$potatoes_requested, collapse = ", "),
      "Potatoes with HMM" = paste(hmm$potatoes_with_genes, collapse = ", "),
      "Genomes" = as.character(hmm$n_genomes),
      "HMM profiles" = as.character(hmm$n_profiles)
    ))
    cli::cli_text("")
    cli::cli_text("Command template:")
    cli::cli_code(hmm$commands$hmmsearch_template)
    cli::cli_text("")
  }

  # Scoring
  if (!is.null(sack@provenance$scoring)) {
    scoring <- sack@provenance$scoring
    cli::cli_h2("Pathway Scoring")
    cli::cli_dl(c(
      "Timestamp" = scoring$timestamp,
      "Pathways scored" = as.character(scoring$n_pathways),
      "Genomes" = as.character(scoring$n_genomes)
    ))
    cli::cli_text("")
    cli::cli_text("Thresholds:")
    cli::cli_ul(c(
      "kofam: {scoring$thresholds$kofam_threshold %||% 'per-gene default'}",
      "blast evalue: {scoring$thresholds$blast_evalue}",
      "blast bitscore: {scoring$thresholds$blast_bitscore}",
      "hmm evalue: {scoring$thresholds$hmm_evalue}"
    ))
  }

  invisible(sack@provenance)
}
