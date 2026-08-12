#' Compute MD5 hash of potato R object
#'
#' Computes a hash of the potato's functional content (nodes, edges, scoring, etc.)
#' to track which version was used for annotation. Excludes metadata fields like
#' json_path and graph cache that don't affect annotation behavior.
#'
#' Specifically excludes cosmetic fields that don't affect detection/scoring:
#' - x, y, x_compounds, y_compounds (visualization coordinates)
#' - notes (documentation only)
#'
#' @param potato Potato S7 object
#' @return Character string with MD5 hash
#' @keywords internal
compute_potato_hash <- function(potato) {
  # Handle both v1 (nodes) and v2 (genes)
  if (S7::S7_inherits(potato, PotatoV2)) {
    # V2 potato
    genes_functional <- lapply(potato@genes, function(gene) {
      # Remove coordinate fields and notes
      gene[!names(gene) %in% c("x", "y", "x_compounds", "y_compounds", "notes")]
    })

    hash_data <- list(
      schema_version = potato@schema_version,
      id = potato@id,
      name = potato@name,
      genes = genes_functional,
      compounds = potato@compounds,
      pathways = potato@pathways,
      tags = potato@tags,
      source = potato@source
    )
  } else {
    # V1 potato
    nodes_functional <- lapply(potato@genes, function(node) {
      # Remove coordinate fields and notes
      node[!names(node) %in% c("x", "y", "x_compounds", "y_compounds", "notes")]
    })

    hash_data <- list(
      id = potato@id,
      name = potato@name,
      nodes = nodes_functional,
      edges = potato@edges,
      tags = potato@tags,
      source = potato@source,
      scoring = potato@scoring,
      input = potato@input,
      output = potato@output
    )
  }

  digest::digest(hash_data, algo = "md5")
}


#' Get potato hashes for multiple potatoes
#'
#' @param potatoes Named list of Potato objects
#' @return Named character vector of hashes
#' @keywords internal
get_potato_hashes <- function(potatoes) {
  hashes <- sapply(potatoes, compute_potato_hash)
  names(hashes) <- names(potatoes)
  hashes
}
