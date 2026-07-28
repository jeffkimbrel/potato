#' Summarize pathway scoring results
#'
#' Provides a high-level summary of pathway scoring results, including
#' detection rates, completeness statistics, and per-pathway summaries.
#'
#' @param sack PotatoSack object with scoring completed
#' @param min_completeness Minimum completeness to consider pathway "present" (default: 1.0)
#' @param verbose Print detailed messages (default: FALSE)
#'
#' @return A list containing:
#'   \item{status}{List with \code{ok} (logical), counts, and detection rate}
#'   \item{summary}{Tibble with per-potato detection statistics}
#'   \item{plot}{ggplot object showing completeness distribution}
#'   \item{messages}{Tibble with \code{type} and \code{message} columns (kable-ready)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sack <- load_potato_sack("results/my_analysis")
#' sack <- annotate_sack(sack)
#' sack <- score_sack(sack)
#'
#' result <- summarize_scoring(sack)
#' result$status
#' result$summary
#' result$plot
#' knitr::kable(result$messages)
#' }
summarize_scoring <- function(sack, min_completeness = 1.0, verbose = FALSE) {

  # Check that scoring has been run
  if (!"scoring" %in% sack@completed_stages) {
    stop("Must run score_sack() before summarizing scoring", call. = FALSE)
  }

  if (is.null(sack@scores)) {
    stop("No scoring results found in sack", call. = FALSE)
  }

  scores <- sack@scores

  # Per-potato summary
  potato_summary <- scores %>%
    dplyr::group_by(potato) %>%
    dplyr::summarise(
      n_genomes = dplyr::n(),
      n_present = sum(present),
      detection_rate = mean(present),
      mean_completeness = mean(completeness),
      median_completeness = stats::median(completeness),
      max_completeness = max(completeness),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(detection_rate))

  # Per-genome summary
  genome_summary <- scores %>%
    dplyr::group_by(genome) %>%
    dplyr::summarise(
      n_potatoes = dplyr::n(),
      n_present = sum(present),
      n_partial = sum(completeness > 0 & completeness < min_completeness),
      mean_completeness = mean(completeness[present]),
      .groups = "drop"
    )

  # Overall completeness distribution
  completeness_dist <- list(
    mean = mean(scores$completeness),
    median = stats::median(scores$completeness),
    q25 = stats::quantile(scores$completeness, 0.25),
    q75 = stats::quantile(scores$completeness, 0.75),
    n_complete = sum(scores$completeness == 1.0),
    n_partial = sum(scores$completeness > 0 & scores$completeness < 1.0),
    n_absent = sum(scores$completeness == 0)
  )

  # Status
  n_present <- sum(scores$present)
  n_total <- nrow(scores)

  status <- list(
    ok = TRUE,
    n_genomes = length(unique(scores$genome)),
    n_potatoes = length(unique(scores$potato)),
    n_present = n_present,
    n_total = n_total,
    detection_rate = n_present / n_total,
    completed_at = if (!is.null(sack@provenance$scoring$timestamp)) {
      sack@provenance$scoring$timestamp
    } else {
      NA
    }
  )

  # Build messages from scoring stage
  sack_msgs <- sack_messages(sack, stage = "scoring", as_dataframe = TRUE)

  if (nrow(sack_msgs) > 0) {
    messages_tbl <- tibble::tibble(
      type = ifelse(sack_msgs$level == "error", "danger",
                   ifelse(sack_msgs$level == "warning", "warning", "info")),
      message = sack_msgs$message
    )
  } else {
    # If no messages collected, create summary messages
    messages_list <- list(
      list(type = "info", message = sprintf("Scored %d genome(s) × %d potato(es)", status$n_genomes, status$n_potatoes)),
      list(type = "success", message = sprintf("Pathways detected: %d / %d (%.1f%%)", n_present, n_total, status$detection_rate * 100))
    )

    if (completeness_dist$n_complete > 0) {
      messages_list[[length(messages_list) + 1]] <- list(
        type = "success",
        message = sprintf("%d pathway(s) with 100%% completeness", completeness_dist$n_complete)
      )
    }

    if (completeness_dist$n_partial > 0) {
      messages_list[[length(messages_list) + 1]] <- list(
        type = "info",
        message = sprintf("%d pathway(s) partially present", completeness_dist$n_partial)
      )
    }

    messages_tbl <- tibble::tibble(
      type = sapply(messages_list, function(x) x$type),
      message = sapply(messages_list, function(x) x$message)
    )
  }

  # Print summary if verbose
  if (verbose) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      cat("\n=== Scoring Summary ===\n")
      cat("Genomes:", status$n_genomes, "\n")
      cat("Potatoes:", status$n_potatoes, "\n")
      cat("Pathways detected:", n_present, "/", n_total,
          sprintf("(%.1f%%)", status$detection_rate * 100), "\n")
      cat("\nCompleteness distribution:\n")
      cat("  Complete (1.0):", completeness_dist$n_complete, "\n")
      cat("  Partial (0-1.0):", completeness_dist$n_partial, "\n")
      cat("  Absent (0.0):", completeness_dist$n_absent, "\n")
    } else {
      cli::cli_h2("Scoring Summary")
      cli::cli_alert_info("Genomes: {status$n_genomes}")
      cli::cli_alert_info("Potatoes: {status$n_potatoes}")
      cli::cli_alert_success("Pathways detected: {n_present} / {n_total} ({round(status$detection_rate * 100, 1)}%)")

      cli::cli_h3("Completeness distribution")
      cli::cli_ul()
      cli::cli_li("Complete (1.0): {completeness_dist$n_complete}")
      cli::cli_li("Partial (0-1.0): {completeness_dist$n_partial}")
      cli::cli_li("Absent (0.0): {completeness_dist$n_absent}")
      cli::cli_end()
    }

    # Top detected potatoes
    if (nrow(potato_summary) > 0) {
      top_n <- min(5, nrow(potato_summary))
      top_potatoes <- potato_summary %>%
        dplyr::slice_head(n = top_n)

      if (!requireNamespace("cli", quietly = TRUE)) {
        cat("\nTop detected potatoes:\n")
        for (i in seq_len(nrow(top_potatoes))) {
          cat("  ", top_potatoes$potato[i], ": ",
              top_potatoes$n_present[i], "/", top_potatoes$n_genomes[i],
              " (", sprintf("%.1f%%", top_potatoes$detection_rate[i] * 100), ")\n", sep = "")
        }
      } else {
        cli::cli_h3("Top detected potatoes")
        for (i in seq_len(nrow(top_potatoes))) {
          cli::cli_text("{cli::symbol$bullet} {top_potatoes$potato[i]}: {top_potatoes$n_present[i]}/{top_potatoes$n_genomes[i]} ({round(top_potatoes$detection_rate[i] * 100, 1)}%)")
        }
      }
    }
  }

  # Create plot: boxplot of completeness per potato, colored by presence
  plot <- NULL

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Determine fill color based on whether any genome has pathway present
    potato_colors <- potato_summary %>%
      dplyr::mutate(
        fill_color = ifelse(n_present > 0, "#2ecc71", "#e74c3c")  # green if present, red if absent
      )

    scores_with_color <- scores %>%
      dplyr::left_join(
        potato_colors %>% dplyr::select(potato, fill_color),
        by = "potato"
      )

    plot <- ggplot2::ggplot(scores_with_color, ggplot2::aes(x = potato, y = completeness, fill = fill_color)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::scale_fill_identity() +
      ggplot2::geom_hline(yintercept = min_completeness, linetype = "dashed", color = "gray40", alpha = 0.5) +
      ggplot2::labs(
        x = "Potato",
        y = "Completeness",
        title = "Pathway completeness across genomes",
        subtitle = "Green = pathway present in at least one genome; Red = absent in all genomes"
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  # Return amplify-style list
  list(
    status = status,
    summary = if (!requireNamespace("tibble", quietly = TRUE)) {
      potato_summary
    } else {
      tibble::as_tibble(potato_summary)
    },
    plot = plot,
    messages = messages_tbl
  )
}


#' Get detailed scoring results
#'
#' Returns the full scoring results table with per-genome, per-potato scores.
#' This provides access to completeness values and detected genes that are
#' used for scoring.
#'
#' @param sack PotatoSack object with scoring completed
#' @param genome Optional filter by genome name
#' @param potato Optional filter by potato ID
#' @param present_only Only return pathways marked as present (default: FALSE)
#'
#' @return Tibble with full scoring results
#' @export
#'
#' @examples
#' \dontrun{
#' sack <- load_potato_sack("results/my_analysis")
#' sack <- score_sack(sack)
#'
#' # Get all scores
#' details <- get_scoring_details(sack)
#'
#' # Get scores for one genome
#' details <- get_scoring_details(sack, genome = "genome001")
#'
#' # Get only present pathways
#' details <- get_scoring_details(sack, present_only = TRUE)
#' }
get_scoring_details <- function(sack, genome = NULL, potato = NULL, present_only = FALSE) {

  if (!"scoring" %in% sack@completed_stages) {
    stop("Must run score_sack() before accessing scoring details", call. = FALSE)
  }

  if (is.null(sack@scores)) {
    stop("No scoring results found in sack", call. = FALSE)
  }

  results <- sack@scores

  # Apply filters
  if (!is.null(genome)) {
    results <- results[results$genome == genome, ]
  }

  if (!is.null(potato)) {
    results <- results[results$potato == potato, ]
  }

  if (present_only) {
    results <- results[results$present, ]
  }

  if (!requireNamespace("tibble", quietly = TRUE)) {
    results
  } else {
    tibble::as_tibble(results)
  }
}


#' Get detected genes for a pathway
#'
#' Returns the list of genes detected for a specific genome-potato combination.
#' Useful for understanding which genes contributed to the pathway score.
#'
#' @param sack PotatoSack object with scoring completed
#' @param genome Genome name
#' @param potato Potato ID
#' @param format Return format: "vector" (default) or "dataframe"
#'
#' @return Character vector or data frame of detected genes
#' @export
#'
#' @examples
#' \dontrun{
#' sack <- load_potato_sack("results/my_analysis")
#' sack <- score_sack(sack)
#'
#' # Get detected genes for one pathway
#' genes <- get_detected_genes(sack, genome = "genome001", potato = "nitrogen_fixation")
#' }
get_detected_genes <- function(sack, genome, potato, format = c("vector", "dataframe")) {

  format <- match.arg(format)

  if (!"scoring" %in% sack@completed_stages) {
    stop("Must run score_sack() before accessing detected genes", call. = FALSE)
  }

  if (is.null(sack@scores)) {
    stop("No scoring results found in sack", call. = FALSE)
  }

  # Find the row
  row <- sack@scores[sack@scores$genome == genome & sack@scores$potato == potato, ]

  if (nrow(row) == 0) {
    stop("No scoring result found for genome '", genome, "' and potato '", potato, "'", call. = FALSE)
  }

  detected <- row$detected_genes[[1]]

  if (format == "vector") {
    return(detected)
  } else {
    # Return as data frame with gene info from potato
    potato_obj <- sack@potatoes[[potato]]
    if (is.null(potato_obj)) {
      warning("Potato '", potato, "' not found in sack")
      return(data.frame(gene_id = detected, stringsAsFactors = FALSE))
    }

    # Build data frame with gene details
    gene_info <- lapply(detected, function(gene_id) {
      node <- potato_obj@nodes[[which(sapply(potato_obj@nodes, function(n) n$id == gene_id))]]
      if (length(node) == 0) {
        data.frame(gene_id = gene_id, name = NA, type = NA, required = NA, stringsAsFactors = FALSE)
      } else {
        data.frame(
          gene_id = gene_id,
          name = if (!is.null(node$name)) node$name else NA,
          type = if (!is.null(node$type)) node$type else NA,
          required = if (!is.null(node$required)) node$required else NA,
          is_marker = if (!is.null(node$is_marker)) node$is_marker else NA,
          stringsAsFactors = FALSE
        )
      }
    })

    result <- do.call(rbind, gene_info)

    if (!requireNamespace("tibble", quietly = TRUE)) {
      result
    } else {
      tibble::as_tibble(result)
    }
  }
}
