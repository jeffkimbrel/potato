#' Summarize genome addition history
#'
#' Displays messages from add_genomes() calls. Messages are reconstructed from
#' provenance data stored in the sack, so this works even if the original
#' add_genomes() call was in an eval=FALSE chunk.
#'
#' @param sack PotatoSack object
#' @param verbose Logical. If TRUE, prints detailed messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{summary}{Tibble with per-call statistics (call_number, timestamp, n_added, total_after)}
#'   \item{messages}{Tibble with columns type, message}
#'   \item{status}{List with ok (logical), total_genomes (integer)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sack <- add_genomes(sack, "~/genomes/*.faa")
#' result <- summarize_add_genomes(sack)
#' result$summary
#' result$messages  # Render this tibble in qmd
#' }
summarize_add_genomes <- function(sack, verbose = FALSE) {

  # Check if genomes have been added
  if (is.null(sack@provenance$genomes) || length(sack@provenance$genomes) == 0) {
    cli::cli_alert_info("No genomes have been added to this sack")
    return(list(
      summary = tibble::tibble(),
      messages = tibble::tibble(type = character(), message = character()),
      status = list(ok = TRUE, total_genomes = 0)
    ))
  }

  # Build summary tibble
  summary_rows <- list()
  cumulative_total <- 0

  for (i in seq_along(sack@provenance$genomes)) {
    prov <- sack@provenance$genomes[[i]]
    cumulative_total <- cumulative_total + prov$n_added

    summary_rows[[i]] <- list(
      call_number = i,
      timestamp = prov$timestamp,
      n_added = prov$n_added,
      total_after = cumulative_total
    )
  }

  summary_tbl <- tibble::tibble(
    call_number = sapply(summary_rows, function(x) x$call_number),
    timestamp = sapply(summary_rows, function(x) x$timestamp),
    n_added = sapply(summary_rows, function(x) x$n_added),
    total_after = sapply(summary_rows, function(x) x$total_after)
  )

  # Reconstruct messages
  messages_list <- list()

  for (i in seq_along(sack@provenance$genomes)) {
    prov <- sack@provenance$genomes[[i]]

    # Info message about what was found
    messages_list[[length(messages_list) + 1]] <- list(
      type = "info",
      message = glue::glue("Call {i} ({prov$timestamp}): Found {prov$n_added} genome file{ifelse(prov$n_added == 1, '', 's')}")
    )

    # Success message
    messages_list[[length(messages_list) + 1]] <- list(
      type = "success",
      message = glue::glue("Added {prov$n_added} genome{ifelse(prov$n_added == 1, '', 's')}. Total: {summary_rows[[i]]$total_after}")
    )
  }

  messages_tbl <- tibble::tibble(
    type = sapply(messages_list, function(x) x$type),
    message = sapply(messages_list, function(x) x$message)
  )

  # Print messages if verbose
  if (verbose) {
    cli::cli_h2("Genome Addition Summary")

    for (i in seq_len(nrow(messages_tbl))) {
      msg <- messages_tbl[i, ]
      switch(msg$type,
        "info" = cli::cli_alert_info(msg$message),
        "success" = cli::cli_alert_success(msg$message),
        "warning" = cli::cli_alert_warning(msg$message),
        "danger" = cli::cli_alert_danger(msg$message),
        cli::cli_text(msg$message)
      )
    }

    cli::cli_text("")
    cli::cli_text("Current sack contains {length(sack@genomes)} genome{?s}")
  }

  # Build status
  status <- list(
    ok = TRUE,
    total_genomes = length(sack@genomes),
    n_calls = length(sack@provenance$genomes)
  )

  # Return summary
  list(
    summary = summary_tbl,
    messages = messages_tbl,
    status = status
  )
}
