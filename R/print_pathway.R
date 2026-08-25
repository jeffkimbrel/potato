#' Print potato structure as text
#'
#' Shows pathway flow with genes, compounds, and structure
#'
#' @param potato Potato S7 object, v2 potato object, or path to JSON
#' @param compact Show compact one-line view (default: TRUE)
#' @param show_compounds Include compound names in flow (default: TRUE, ignored if compact)
#' @param show_databases Show database annotations (kofam, blast, hmm, etc.) (default: FALSE)
#' @param show_ec Show EC numbers (default: FALSE)
#'
#' @export

print_potato <- function(potato, compact = TRUE, show_compounds = TRUE, show_databases = FALSE, show_ec = FALSE) {

  # Load if path provided
  if (is.character(potato)) {
    potato <- load_potato_v2(potato)
  }

  cli::cli_h1(potato@name)
  cli::cli_text("Source: {.field {potato@source}}")
  cli::cli_text("Tags: {.val {paste(potato@tags, collapse = ', ')}}")
  if (!is.null(potato@notes) && nchar(potato@notes) > 0) {
    cli::cli_text("Notes: {.field {potato@notes}}")
  }
  cli::cli_rule()

  # Show genes
  cli::cli_h2("Genes ({length(potato@genes)})")
  for (gene in potato@genes) {
    cli::cli_text("{gene$name} [{.field {gene$id}}]")

    if (show_ec && !is.null(gene$ec) && length(gene$ec) > 0) {
      cli::cli_text("  {.dim EC: {paste(gene$ec, collapse = ', ')}}")
    }

    if (show_databases && !is.null(gene$databases)) {
      db_info <- character()
      if (!is.null(gene$databases$kofam)) {
        ko_terms <- unlist(gene$databases$kofam)
        db_info <- c(db_info, paste0("KOfam: ", paste(ko_terms, collapse = ", ")))
      }
      if (!is.null(gene$databases$blast)) {
        blast_terms <- unlist(gene$databases$blast)
        db_info <- c(db_info, paste0("BLAST: ", paste(blast_terms, collapse = ", ")))
      }
      if (!is.null(gene$databases$hmm)) {
        hmm_terms <- unlist(gene$databases$hmm)
        db_info <- c(db_info, paste0("HMM: ", paste(hmm_terms, collapse = ", ")))
      }

      if (length(db_info) > 0) {
        cli::cli_text("  {.dim ({paste(db_info, collapse = ', ')})}")
      }
    }

    if (!is.null(gene$reactions) && length(gene$reactions) > 0) {
      cli::cli_text("  {.dim Reactions: {paste(gene$reactions, collapse = ', ')}}")
    }
  }

  cli::cli_text("")

  # Show compounds
  if (!is.null(potato@compounds) && length(potato@compounds) > 0) {
    cli::cli_h2("Compounds ({length(potato@compounds)})")
    for (compound in potato@compounds) {
      cli::cli_text("{compound$name} [{.field {compound$id}}]")
    }
    cli::cli_text("")
  }

  # Show pathway(s)
  if (!is.null(potato@pathways)) {
    pathway_names <- names(potato@pathways)
    cli::cli_h2("Pathway{?s} ({length(pathway_names)})")

    for (pathway_id in pathway_names) {
      pathway <- potato@pathways[[pathway_id]]
      cli::cli_text("{.strong {pathway$name %||% pathway_id}}")

      if (!is.null(pathway$notes) && nchar(pathway$notes) > 0) {
        cli::cli_text("  {.dim {pathway$notes}}")
      }

      # Show edges
      if (!is.null(pathway$edges) && length(pathway$edges) > 0) {
        cli::cli_text("  {.dim {length(pathway$edges)} edges}")
        for (edge in head(pathway$edges, 3)) {
          from_label <- if (grepl("^C\\d+", edge$from)) {
            compound <- Find(function(c) c$id == edge$from, potato@compounds)
            if (!is.null(compound)) compound$name else edge$from
          } else {
            edge$from
          }
          to_label <- if (grepl("^C\\d+", edge$to)) {
            compound <- Find(function(c) c$id == edge$to, potato@compounds)
            if (!is.null(compound)) compound$name else edge$to
          } else {
            edge$to
          }
          rxn_str <- if (!is.null(edge$reaction)) {
            rxn_display <- if (length(edge$reaction) > 1) paste(edge$reaction, collapse = ", ") else edge$reaction
            paste0(" [", rxn_display, "]")
          } else ""
          cli::cli_text("    {from_label} → {to_label}{rxn_str}")
        }
        if (length(pathway$edges) > 3) {
          cli::cli_text("    {.dim ... and {length(pathway$edges) - 3} more}")
        }
      }

      cli::cli_text("")
    }
  }

  return(invisible(potato))
}
