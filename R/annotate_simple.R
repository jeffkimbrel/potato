#' Simple annotation - just produce gene-level results
#'
#' @param sack PotatoSack object
#' @export
annotate_sack_simple <- function(sack) {

  # Check jakomics available
  if (!exists("jakomics")) {
    stop("jakomics not loaded", call. = FALSE)
  }

  # Get BLAST config
  if (!"blast" %in% names(sack@config$databases)) {
    stop("No BLAST database configured", call. = FALSE)
  }

  blast_config <- sack@config$databases$blast
  blast_files <- blast_config$files %||% blast_config$path
  blast_db <- blast_files[1]  # Use first file for now

  # Build term → potato gene mapping
  # term_map[[term]] = list(potato_id, gene_id, gene_name)
  term_map <- list()

  for (potato in sack@potatoes) {
    for (gene in potato@genes) {
      if (is.null(gene$type) || gene$type != "enzyme") next

      # Get BLAST terms for this gene
      if (!is.null(gene$databases) && "blast" %in% names(gene$databases)) {
        terms <- gene$databases$blast
        for (term in terms) {
          term_map[[term]] <- list(
            potato_id = potato@id,
            gene_id = gene$id,
            gene_name = gene$name %||% gene$id
          )
        }
      }
    }
  }

  if (length(term_map) == 0) {
    message("No BLAST terms found in potatoes")
    return(data.frame())
  }

  message("Found ", length(term_map), " BLAST terms across ", length(sack@potatoes), " potato(es)")

  # Run BLAST on each genome - return nested tibble
  genome_results <- list()

  for (genome in sack@genomes) {
    message("Running BLAST: ", genome$short_name)

    # Call jakomics directly
    blast_results <- jakomics$blast$run_blast(
      type = "prot",
      q = genome$file_path,
      db = blast_db,
      e = 1e-7,
      make = FALSE,
      return_query_results = FALSE,  # Return by subject (database term)
      echo = FALSE
    )

    # Parse results for this genome
    gene_hits <- list()

    for (term in names(blast_results)) {
      # Check if this term is in our potato mapping
      if (!term %in% names(term_map)) next

      hits <- blast_results[[term]]
      potato_info <- term_map[[term]]

      # Each hit is a match
      for (hit in hits) {
        gene_hits <- c(gene_hits, list(data.frame(
          potato_id = potato_info$potato_id,
          gene = potato_info$gene_id,
          product = potato_info$gene_name,
          type = "BLAST",
          id = term,
          locus_tag = hit$query,
          score = hit$bit_score,
          evalue = hit$eval,
          pident = hit$percent,
          stringsAsFactors = FALSE
        )))
      }
    }

    # Convert to tibble for this genome
    if (length(gene_hits) > 0) {
      genome_tibble <- dplyr::bind_rows(gene_hits)
    } else {
      genome_tibble <- tibble::tibble(
        potato_id = character(),
        gene = character(),
        product = character(),
        type = character(),
        id = character(),
        locus_tag = character(),
        score = numeric(),
        evalue = numeric(),
        pident = numeric()
      )
    }

    genome_results[[genome$short_name]] <- genome_tibble
  }

  # Create nested tibble
  result <- tibble::tibble(
    genome = names(genome_results),
    genes = genome_results
  )

  total_hits <- sum(sapply(genome_results, nrow))
  message("Found ", total_hits, " gene matches across ", nrow(result), " genome(s)")

  result
}
