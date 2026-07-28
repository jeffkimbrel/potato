#' Map annotation hits to potato genes
#'
#' Takes raw annotation results and maps them to specific genes in potatoes.
#' Adds columns identifying which potato(s) and gene(s) each hit corresponds to.
#'
#' @param hits Data frame of annotation results (from unnest_all_annotations())
#' @param potatoes List of Potato objects or single Potato
#' @param tool Character. Tool type ("kofam", "blast", "pfam", "hmm")
#'
#' @returns Data frame with additional columns:
#'   - potato_id: which potato contains this term
#'   - potato_gene: which gene in the potato
#'   - potato_gene_name: gene name from potato
#'   - is_marker: logical, is this a marker gene
#'
#' @details
#' A single hit can map to multiple potatoes if the same term appears in
#' multiple pathways. The result will have one row per hit-potato combination.
#'
#' @keywords internal
map_hits_to_potatoes <- function(hits, potatoes, tool = "kofam") {

  if (inherits(potatoes, "Potato")) {
    potatoes <- list(potatoes)
  }

  if (nrow(hits) == 0) {
    return(hits)
  }

  # Build lookup: term -> list of (potato_id, gene_id, gene_name, is_marker)
  term_col <- switch(tool,
    kofam = "ko",
    pfam = "pfam",
    blast = "subject",  # BLAST results use 'subject' for database gene ID
    hmm = "target",     # HMM results typically use 'target'
    stop("Unknown tool type: ", tool)
  )

  if (!term_col %in% names(hits)) {
    warning("No '", term_col, "' column in hits data frame for tool '", tool, "'")
    return(hits)
  }

  # Build mapping from annotation term -> potato genes
  mapping_list <- list()

  for (potato in potatoes) {
    for (node in potato@nodes) {
      if (is.null(node$type) || node$type != "enzyme") next

      terms <- character(0)

      # New schema: extract terms from databases field
      if (!is.null(node$databases)) {
        for (db_name in names(node$databases)) {
          # For now, accept any database (we'll refine this later to match tool type)
          db_terms <- node$databases[[db_name]]
          terms <- c(terms, if (is.list(db_terms)) unlist(db_terms) else db_terms)
        }
      }

      # Legacy schema: Get detection terms for this tool
      legacy_terms <- switch(tool,
        kofam = node$ko,
        pfam = node$pfam,
        blast = node$blast_terms,
        hmm = node$hmm,
        NULL
      )

      if (!is.null(legacy_terms)) {
        terms <- c(terms, if (is.list(legacy_terms)) unlist(legacy_terms) else legacy_terms)
      }

      if (length(terms) == 0) next

      # Each term maps to this gene
      for (term in terms) {
        if (!term %in% names(mapping_list)) {
          mapping_list[[term]] <- list()
        }

        mapping_list[[term]] <- c(mapping_list[[term]], list(list(
          potato_id = potato@id,
          potato_gene = node$id,
          potato_gene_name = if (!is.null(node$name)) node$name else node$id,
          is_marker = if (!is.null(node$marker)) node$marker else FALSE
        )))
      }
    }
  }

  # Map each hit
  result_rows <- list()

  for (i in seq_len(nrow(hits))) {
    hit <- hits[i, ]
    term <- hit[[term_col]]

    # Look up this term
    if (term %in% names(mapping_list)) {
      # One row per potato that contains this term
      for (gene_info in mapping_list[[term]]) {
        new_row <- hit
        new_row$potato_id <- gene_info$potato_id
        new_row$potato_gene <- gene_info$potato_gene
        new_row$potato_gene_name <- gene_info$potato_gene_name
        new_row$is_marker <- gene_info$is_marker
        result_rows <- c(result_rows, list(new_row))
      }
    } else {
      # Term not in any potato (shouldn't happen if we searched correctly)
      new_row <- hit
      new_row$potato_id <- NA_character_
      new_row$potato_gene <- NA_character_
      new_row$potato_gene_name <- NA_character_
      new_row$is_marker <- FALSE
      result_rows <- c(result_rows, list(new_row))
    }
  }

  do.call(rbind, result_rows)
}


#' Map annotation table to potatoes
#'
#' Maps all results in an annotation table to potato genes
#'
#' @param anno_table Annotation table from initialize_annotation_table()
#' @param potatoes List of Potato objects
#' @param config Potato config object (to look up database types)
#' @param databases Character vector of which databases to map. NULL = all
#'
#' @returns Modified annotation table with mapped results
#' @keywords internal
map_annotation_table <- function(anno_table, potatoes, config, databases = NULL) {

  if (!inherits(anno_table, "potato_annotation_table")) {
    stop("Input must be a potato_annotation_table", call. = FALSE)
  }

  db_cols <- setdiff(names(anno_table), c("genome", "file_object"))

  if (!is.null(databases)) {
    db_cols <- intersect(db_cols, databases)
  }

  # Map each database column
  for (db_name in db_cols) {
    # Look up database type from config
    if (!db_name %in% names(config$databases)) {
      warning("Database '", db_name, "' not found in config, skipping mapping")
      next
    }

    db_type <- config$databases[[db_name]]$type

    for (i in seq_len(nrow(anno_table))) {
      result <- anno_table[[db_name]][[i]]

      if (is.null(result) || !is.data.frame(result) || nrow(result) == 0) {
        next
      }

      # Map this genome's results using the database type
      mapped <- map_hits_to_potatoes(result, potatoes, tool = db_type)
      anno_table[[db_name]][[i]] <- mapped
    }
  }

  anno_table
}
