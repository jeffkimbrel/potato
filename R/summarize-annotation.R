#' Summarize annotation results
#'
#' Provides a high-level summary of annotation stage results, including
#' database statistics, tool performance, and any warnings or errors.
#'
#' @param sack PotatoSack object with annotation completed
#' @param verbose Print detailed messages (default: FALSE)
#'
#' @return A list containing:
#'   \item{status}{List with \code{ok} (logical), summary stats, and timestamps}
#'   \item{summary}{Tibble with per-genome annotation statistics}
#'   \item{plot}{ggplot object showing hit distributions (NULL for now)}
#'   \item{messages}{Tibble with \code{type} and \code{message} columns (kable-ready)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sack <- load_potato_sack("results/my_analysis")
#' sack <- annotate_sack(sack)
#' result <- summarize_annotation(sack)
#'
#' result$status$ok
#' result$summary
#' knitr::kable(result$messages)
#' }
summarize_annotation <- function(sack, verbose = FALSE) {

  # Check that annotation has been run
  if (!"annotation" %in% sack@completed_stages) {
    stop("Must run annotate_sack() before summarizing annotation", call. = FALSE)
  }

  if (is.null(sack@results)) {
    stop("No annotation results found in sack", call. = FALSE)
  }

  results_df <- sack@results
  db_cols <- setdiff(names(results_df), c("genome", "file_object"))

  # Per-genome summary
  genome_summary <- data.frame(
    genome = results_df$genome,
    stringsAsFactors = FALSE
  )

  # Count hits per database
  for (db in db_cols) {
    col_name <- paste0("n_", db)
    genome_summary[[col_name]] <- sapply(results_df[[db]], function(x) {
      if (is.null(x) || nrow(x) == 0) 0 else nrow(x)
    })
  }

  # Total hits
  genome_summary$n_total_hits <- rowSums(genome_summary[, paste0("n_", db_cols), drop = FALSE])

  # Per-database summary for display
  db_names <- names(sack@config$databases)
  db_types <- sapply(sack@config$databases, function(db) db$type)

  db_summary <- data.frame(
    database = db_names,
    type = db_types,
    n_hits = 0,
    stringsAsFactors = FALSE
  )

  # Add hit counts per database
  for (i in seq_len(nrow(db_summary))) {
    db_name <- db_summary$database[i]
    if (db_name %in% names(results_df)) {
      db_summary$n_hits[i] <- sum(sapply(results_df[[db_name]], function(x) {
        if (is.null(x) || nrow(x) == 0) 0 else nrow(x)
      }))
    }
  }

  # Get collected messages from sack
  sack_msgs <- sack_messages(sack, as_dataframe = TRUE)
  anno_msgs <- if (nrow(sack_msgs) > 0) {
    sack_msgs[sack_msgs$stage %in% c("kofam", "blast", "hmm", "pfam", "annotation"), ]
  } else {
    data.frame(level = character(), message = character(), stringsAsFactors = FALSE)
  }

  # Convert to amplify-style messages (type, message)
  if (nrow(anno_msgs) > 0) {
    messages_tbl <- tibble::tibble(
      type = ifelse(anno_msgs$level == "error", "danger",
                   ifelse(anno_msgs$level == "warning", "warning", "info")),
      message = anno_msgs$message
    )
  } else {
    messages_tbl <- tibble::tibble(
      type = character(),
      message = character()
    )
  }

  # Status
  n_errors <- sum(anno_msgs$level == "error")
  n_warnings <- sum(anno_msgs$level == "warning")

  status <- list(
    ok = n_errors == 0,
    n_genomes = nrow(results_df),
    n_databases = length(db_names),
    n_warnings = n_warnings,
    n_errors = n_errors,
    total_hits = sum(genome_summary$n_total_hits),
    completed_at = if (!is.null(sack@provenance$annotation$timestamp)) {
      sack@provenance$annotation$timestamp
    } else {
      NA
    }
  )

  # Print summary if verbose
  if (verbose) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      cat("\n=== Annotation Summary ===\n")
      cat("Genomes:", status$n_genomes, "\n")
      cat("Databases:", status$n_databases, "\n")
      cat("Total hits:", status$total_hits, "\n")
      cat("Status:", if (status$ok) "OK" else "ERRORS", "\n")
      if (n_warnings > 0) cat("Warnings:", n_warnings, "\n")
      if (n_errors > 0) cat("Errors:", n_errors, "\n")
    } else {
      cli::cli_h2("Annotation Summary")
      cli::cli_alert_info("Genomes: {status$n_genomes}")
      cli::cli_alert_info("Databases: {status$n_databases}")
      cli::cli_alert_info("Total hits: {status$total_hits}")

      if (status$ok && n_warnings == 0) {
        cli::cli_alert_success("Status: OK")
      } else if (n_errors > 0) {
        cli::cli_alert_danger("Status: {n_errors} error{?s}")
      } else {
        cli::cli_alert_warning("Status: {n_warnings} warning{?s}")
      }
    }

    # Show database summary
    if (!requireNamespace("cli", quietly = TRUE)) {
      cat("\nDatabase hits:\n")
      for (i in seq_len(nrow(db_summary))) {
        cat("  ", db_summary$database[i], " (", db_summary$type[i], "): ",
            db_summary$n_hits[i], "\n", sep = "")
      }
    } else {
      cli::cli_h3("Database hits")
      for (i in seq_len(nrow(db_summary))) {
        cli::cli_text("{cli::symbol$bullet} {db_summary$database[i]} ({db_summary$type[i]}): {db_summary$n_hits[i]} hit{?s}")
      }
    }
  }

  # Create plot: boxplot of hits per potato by database, colored by passed status
  plot <- NULL

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    all_hits <- get_annotation_details(sack)

    if (nrow(all_hits) > 0 && "potato_id" %in% names(all_hits) && "database" %in% names(all_hits)) {
      # Select unique hits and count by genome, database, potato, passed status
      plot_data <- all_hits %>%
        dplyr::select(genome, potato_id, potato_gene, database, passed) %>%
        dplyr::distinct() %>%
        dplyr::count(genome, database, potato_id, passed)

      if (nrow(plot_data) > 0) {
        plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = database, y = n, color = passed)) +
          ggplot2::geom_boxplot() +
          ggplot2::facet_wrap(~potato_id) +
          ggplot2::scale_color_manual(
            values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c"),
            labels = c("TRUE" = "Passed", "FALSE" = "Failed"),
            name = "Threshold"
          ) +
          ggplot2::labs(
            x = "Database",
            y = "Number of hits per genome",
            title = "Annotation hits by potato and threshold status"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }
    }
  }

  # Return amplify-style list
  list(
    status = status,
    summary = if (!requireNamespace("tibble", quietly = TRUE)) {
      genome_summary
    } else {
      tibble::as_tibble(genome_summary)
    },
    plot = plot,
    messages = messages_tbl
  )
}


#' Get detailed annotation results
#'
#' Returns the full annotation results table with all hits from all tools.
#' This provides access to per-hit details that are not included in the
#' simplified summary.
#'
#' @param sack PotatoSack object with annotation completed
#' @param genome Optional filter by genome name
#' @param potato Optional filter by potato ID
#' @param tool Optional filter by tool type (e.g., "kofam", "blast")
#'
#' @return Tibble with full annotation results
#' @export
#'
#' @examples
#' \dontrun{
#' sack <- load_potato_sack("results/my_analysis")
#' sack <- annotate_sack(sack)
#'
#' # Get all results
#' details <- get_annotation_details(sack)
#'
#' # Get results for one genome
#' details <- get_annotation_details(sack, genome = "genome001")
#'
#' # Get kofam results only
#' details <- get_annotation_details(sack, tool = "kofam")
#' }
get_annotation_details <- function(sack, genome = NULL, potato = NULL, tool = NULL) {

  if (!"annotation" %in% sack@completed_stages) {
    stop("Must run annotate_sack() before accessing annotation details", call. = FALSE)
  }

  if (is.null(sack@results)) {
    stop("No annotation results found in sack", call. = FALSE)
  }

  results <- sack@results
  db_cols <- setdiff(names(results), c("genome", "file_object"))

  # Filter by genome first
  if (!is.null(genome)) {
    results <- results[results$genome == genome, ]
  }

  # If tool specified, filter to databases of that type
  if (!is.null(tool)) {
    # Find databases matching this tool type
    matching_dbs <- names(sack@config$databases)[
      sapply(sack@config$databases, function(db) db$type == tool)
    ]
    db_cols <- intersect(db_cols, matching_dbs)

    if (length(db_cols) == 0) {
      warning("No databases found for tool type '", tool, "'")
      return(data.frame())
    }
  }

  # Unnest each database and combine
  all_hits <- list()
  for (db in db_cols) {
    for (i in seq_len(nrow(results))) {
      hits <- results[[db]][[i]]
      if (!is.null(hits) && nrow(hits) > 0) {
        # Database column should already be in hits from annotation
        # If not, add it
        if (!"database" %in% names(hits)) {
          hits$database <- db
        }
        all_hits <- c(all_hits, list(hits))
      }
    }
  }

  if (length(all_hits) == 0) {
    return(data.frame())
  }

  combined <- dplyr::bind_rows(all_hits)

  # Filter by potato if specified (only if potato column exists)
  if (!is.null(potato) && "potato" %in% names(combined)) {
    combined <- combined[combined$potato == potato, ]
  }

  if (!requireNamespace("tibble", quietly = TRUE)) {
    combined
  } else {
    tibble::as_tibble(combined)
  }
}
