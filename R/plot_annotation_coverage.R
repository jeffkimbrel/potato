#' Plot annotation coverage heatmap
#'
#' Shows which pathways have been annotated with which databases (kofam, blast, hmm).
#' Helps identify gaps in annotation coverage.
#'
#' Cell colors:
#' - Green: Pathway was checked with this database (pathway has genes AND was in potatoes_requested)
#' - Gray: Pathway has no genes for this database (not applicable)
#' - Red: Pathway HAS genes for this database but was NOT checked (annotation gap!)
#'
#' @param sack PotatoSack object with provenance data
#'
#' @returns ggplot2 heatmap
#' @export
plot_annotation_coverage <- function(sack) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} required for plotting")
  }

  # Check if any annotation has been run
  if (is.null(sack@provenance$kofam) &&
      is.null(sack@provenance$blast) &&
      is.null(sack@provenance$hmm)) {
    cli::cli_abort("No annotation provenance found. Run annotation tools first.")
  }

  # Build coverage matrix
  rows <- list()

  for (potato in sack@potatoes) {
    # Determine if this is a multi-pathway network (> 1 pathway)
    is_multi_pathway <- !is.null(potato@pathways) &&
                        is.list(potato@pathways) &&
                        length(potato@pathways) > 1

    if (is_multi_pathway) {
      # Multi-pathway network: one row per pathway
      for (pathway_id in names(potato@pathways)) {
        pathway <- potato@pathways[[pathway_id]]
        pathway_label <- paste0(potato@id, ":", pathway_id)

        # Check which databases this pathway's genes use
        # Extract gene IDs from edges (filter out compounds which start with C)
        pathway_genes <- unique(c(
          sapply(pathway$edges, function(e) e$from),
          sapply(pathway$edges, function(e) e$to)
        ))
        # Remove compounds (start with C followed by digits) and NULLs
        pathway_genes <- pathway_genes[!grepl("^C\\d+", pathway_genes)]
        pathway_genes <- pathway_genes[!is.na(pathway_genes) & nchar(pathway_genes) > 0]

        has_kofam <- FALSE
        has_blast <- FALSE
        has_hmm <- FALSE

        for (gene_id in pathway_genes) {
          # Find gene by ID (genes is a list, not named)
          gene <- Find(function(g) g$id == gene_id, potato@genes)
          if (is.null(gene)) next

          if (!is.null(gene$databases$kofam) && length(gene$databases$kofam) > 0) {
            has_kofam <- TRUE
          }
          if (!is.null(gene$databases$blast) && length(gene$databases$blast) > 0) {
            has_blast <- TRUE
          }
          if (!is.null(gene$databases$hmm) && length(gene$databases$hmm) > 0) {
            has_hmm <- TRUE
          }
        }

        # Check if this potato was checked with each database
        checked_kofam <- !is.null(sack@provenance$kofam) &&
                        potato@id %in% sack@provenance$kofam$potatoes_requested
        checked_blast <- !is.null(sack@provenance$blast) &&
                        potato@id %in% sack@provenance$blast$potatoes_requested
        checked_hmm <- !is.null(sack@provenance$hmm) &&
                      potato@id %in% sack@provenance$hmm$potatoes_requested

        # Determine colors for each database
        kofam_status <- if (!has_kofam) "none" else if (checked_kofam) "checked" else "gap"
        blast_status <- if (!has_blast) "none" else if (checked_blast) "checked" else "gap"
        hmm_status <- if (!has_hmm) "none" else if (checked_hmm) "checked" else "gap"

        rows[[length(rows) + 1]] <- list(
          pathway = pathway_label,
          pathway_name = pathway$name %||% pathway_id,
          database = "kofam",
          status = kofam_status
        )
        rows[[length(rows) + 1]] <- list(
          pathway = pathway_label,
          pathway_name = pathway$name %||% pathway_id,
          database = "blast",
          status = blast_status
        )
        rows[[length(rows) + 1]] <- list(
          pathway = pathway_label,
          pathway_name = pathway$name %||% pathway_id,
          database = "hmm",
          status = hmm_status
        )
      }
    } else {
      # Single pathway potato
      pathway_label <- potato@id

      # Check which databases this potato's genes use
      has_kofam <- FALSE
      has_blast <- FALSE
      has_hmm <- FALSE

      for (gene in potato@genes) {
        if (!is.null(gene$databases$kofam) && length(gene$databases$kofam) > 0) {
          has_kofam <- TRUE
        }
        if (!is.null(gene$databases$blast) && length(gene$databases$blast) > 0) {
          has_blast <- TRUE
        }
        if (!is.null(gene$databases$hmm) && length(gene$databases$hmm) > 0) {
          has_hmm <- TRUE
        }
      }

      # Check if this potato was checked with each database
      checked_kofam <- !is.null(sack@provenance$kofam) &&
                      potato@id %in% sack@provenance$kofam$potatoes_requested
      checked_blast <- !is.null(sack@provenance$blast) &&
                      potato@id %in% sack@provenance$blast$potatoes_requested
      checked_hmm <- !is.null(sack@provenance$hmm) &&
                    potato@id %in% sack@provenance$hmm$potatoes_requested

      # Determine colors for each database
      kofam_status <- if (!has_kofam) "none" else if (checked_kofam) "checked" else "gap"
      blast_status <- if (!has_blast) "none" else if (checked_blast) "checked" else "gap"
      hmm_status <- if (!has_hmm) "none" else if (checked_hmm) "checked" else "gap"

      rows[[length(rows) + 1]] <- list(
        pathway = pathway_label,
        pathway_name = potato@name,
        database = "kofam",
        status = kofam_status
      )
      rows[[length(rows) + 1]] <- list(
        pathway = pathway_label,
        pathway_name = potato@name,
        database = "blast",
        status = blast_status
      )
      rows[[length(rows) + 1]] <- list(
        pathway = pathway_label,
        pathway_name = potato@name,
        database = "hmm",
        status = hmm_status
      )
    }
  }

  # Convert to tibble
  coverage_df <- dplyr::bind_rows(rows)

  # Create plot
  p <- ggplot2::ggplot(coverage_df, ggplot2::aes(x = database, y = pathway, fill = status)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c(
        "checked" = "#2165de",  #
        "none" = "gray40",     # Gray
        "gap" = "#E74C3C"       # Red
      ),
      labels = c(
        "checked" = "Checked",
        "none" = "No genes",
        "gap" = "Gap (not checked)"
      ),
      name = "Status"
    ) +
    ggplot2::labs(
      title = "Annotation Coverage by Database",
      x = "Database",
      y = "Pathway"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(size = 10, angle = 0),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  p
}
