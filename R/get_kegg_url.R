#' Generate KEGG URL for a specific pathway
#'
#' Creates a KEGG database URL with all unique KO IDs and compound IDs
#' from a specific pathway. Useful for viewing pathway entries together.
#'
#' @param potato Potato object or path to potato JSON
#' @param pathway_id Pathway ID within the potato (required for multi-pathway networks)
#' @param genes Logical. Include KO IDs from genes (default: TRUE)
#' @param compounds Logical. Include compound IDs (default: TRUE)
#'
#' @return Character string with KEGG URL, or NULL if no IDs found
#' @export
get_kegg_url <- function(potato, pathway_id = NULL, genes = TRUE, compounds = TRUE) {

  # Load if path provided
  if (is.character(potato)) {
    potato <- load_potato_v2(potato)
  }

  # Handle pathway selection
  if (is.null(pathway_id)) {
    if (length(potato@pathways) > 1) {
      cli::cli_abort(c(
        "Multi-pathway network requires pathway_id parameter",
        "i" = "Available pathways: {paste(names(potato@pathways), collapse=', ')}"
      ))
    }
    pathway_id <- names(potato@pathways)[1]
  }

  if (!pathway_id %in% names(potato@pathways)) {
    cli::cli_abort("Pathway '{pathway_id}' not found")
  }

  pathway <- potato@pathways[[pathway_id]]

  # Collect IDs
  ids <- character(0)

  if (genes) {
    # Get genes used in this pathway
    pathway_gene_ids <- character()

    # From explicit genes array
    if (!is.null(pathway$genes)) {
      pathway_gene_ids <- pathway$genes
    }

    # From edges (gene nodes)
    compound_ids <- sapply(potato@compounds, function(c) c$id)
    for (edge in pathway$edges) {
      if (!edge$from %in% compound_ids) {
        pathway_gene_ids <- c(pathway_gene_ids, edge$from)
      }
      if (!edge$to %in% compound_ids) {
        pathway_gene_ids <- c(pathway_gene_ids, edge$to)
      }
    }
    pathway_gene_ids <- unique(pathway_gene_ids)

    # Extract KO IDs from these genes
    for (gene_id in pathway_gene_ids) {
      gene <- Find(function(g) g$id == gene_id, potato@genes)
      if (!is.null(gene) && !is.null(gene$databases) && !is.null(gene$databases$kofam)) {
        ids <- c(ids, unlist(gene$databases$kofam))
      }
    }
  }

  if (compounds) {
    # Get compounds from pathway edges
    compound_ids_in_pathway <- character()
    for (edge in pathway$edges) {
      compound_ids_in_pathway <- c(compound_ids_in_pathway, edge$from, edge$to)
    }
    compound_ids_in_pathway <- unique(compound_ids_in_pathway)

    # Extract KEGG compound IDs (strip decorations like _internal, _external)
    for (cid in compound_ids_in_pathway) {
      # Extract C##### from decorated IDs like C00160_internal
      if (grepl("^C\\d+", cid)) {
        base_id <- sub("^(C\\d+).*", "\\1", cid)
        ids <- c(ids, base_id)
      }
    }
  }

  # Get unique IDs and sort
  ids <- unique(ids)
  ids <- sort(ids)

  # Build URL
  if (length(ids) == 0) {
    cli::cli_alert_warning("No KEGG IDs found in sack")
    return(NULL)
  }

  url <- paste0("https://www.genome.jp/entry/", paste(ids, collapse = "+"))

  cli::cli_alert_info("Generated KEGG URL with {length(ids)} IDs ({sum(grepl('^K', ids))} KOs, {sum(grepl('^C', ids))} compounds)")

  url
}
