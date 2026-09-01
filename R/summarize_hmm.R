#' Summarize hmm annotation
#'
#' Displays messages from run_hmm(). Messages are retrieved from
#' provenance data stored in the sack, so this works even if the original
#' run_hmm() call was in an eval=FALSE chunk.
#'
#' @param sack PotatoSack object
#' @param verbose Logical. If TRUE, prints detailed messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{summary}{Tibble with annotation summary (n_genomes, n_potatoes, n_profiles, tool_version)}
#'   \item{messages}{Tibble with columns type, message}
#'   \item{status}{List with ok (logical)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sack <- run_hmm(sack)
#' result <- summarize_hmm(sack)
#' result$summary
#' result$messages  # Render this tibble in qmd
#' }
summarize_hmm <- function(sack, verbose = FALSE) {

  # Check if hmm has been run
  if (is.null(sack@provenance$hmm)) {
    cli::cli_alert_info("HMM annotation has not been run on this sack")
    return(list(
      summary = tibble::tibble(),
      messages = tibble::tibble(type = character(), message = character()),
      status = list(ok = FALSE)
    ))
  }

  prov <- sack@provenance$hmm

  # Build summary tibble
  summary_tbl <- tibble::tibble(
    timestamp = prov$timestamp,
    n_genomes = prov$n_genomes,
    n_potatoes = length(prov$potatoes_requested),
    n_profiles = prov$n_profiles,
    tool_version = prov$tool_version,
    conda_env = prov$conda_env %||% "none",
    workers = prov$workers
  )

  # Get messages from provenance (may not exist if run_hmm was called before message tracking)
  if (!is.null(prov$messages) && nrow(prov$messages) > 0) {
    messages_tbl <- prov$messages
  } else {
    # Reconstruct basic messages from provenance
    messages_tbl <- tibble::tibble(
      type = c("info", "success"),
      message = c(
        glue::glue("HMM annotation run on {prov$n_genomes} genome{ifelse(prov$n_genomes == 1, '', 's')}"),
        glue::glue("Annotation complete: {prov$n_profiles} HMM profiles across {prov$n_potatoes} potato{ifelse(prov$n_potatoes == 1, '', 'es')}")
      )
    )
  }

  # Print messages if verbose
  if (verbose) {
    cli::cli_h2("HMM Annotation Summary")

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
    cli::cli_dl(c(
      "Genomes" = as.character(prov$n_genomes),
      "Potatoes" = paste(prov$potatoes_requested, collapse = ", "),
      "HMM profiles" = as.character(prov$n_profiles),
      "Tool version" = prov$tool_version
    ))
  }

  # Build status
  status <- list(
    ok = TRUE
  )

  # Return summary
  list(
    summary = summary_tbl,
    messages = messages_tbl,
    status = status
  )
}
