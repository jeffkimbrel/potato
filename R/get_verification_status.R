#' Get verification status of pathways
#'
#' Shows which pathways in potatoes are verified or need verification.
#' Can check a single potato or all potatoes in a sack.
#'
#' @param potato_or_sack Potato object, potato path, or PotatoSack object
#'
#' @return Prints verification summary and invisibly returns a tibble
#' @export
get_verification_status <- function(potato_or_sack) {

  # Determine if sack or single potato
  if (inherits(potato_or_sack, "potato::PotatoSack")) {
    potatoes <- S7::prop(potato_or_sack, "potatoes")
  } else {
    # Single potato (object or path)
    if (is.character(potato_or_sack)) {
      potato <- load_potato(potato_or_sack)
    } else {
      potato <- potato_or_sack
    }
    potatoes <- list(potato)
    names(potatoes) <- potato@id
  }

  # Collect results
  results <- list()

  for (potato_id in names(potatoes)) {
    potato <- potatoes[[potato_id]]

    # Check if multi-pathway
    is_network <- !is.null(potato@edges) &&
                  is.list(potato@edges) &&
                  length(names(potato@edges)) > 0 &&
                  !is.null(potato@edges[[1]]$type)

    if (is_network) {
      # Multi-pathway: check each pathway
      pathway_status <- sapply(names(potato@edges), function(pw_id) {
        potato@edges[[pw_id]]$verified %||% FALSE
      })

      n_verified <- sum(pathway_status)
      n_total <- length(pathway_status)
      unverified_names <- names(pathway_status)[!pathway_status]

      results[[potato_id]] <- list(
        verified = n_verified,
        total = n_total,
        unverified = if (length(unverified_names) > 0) unverified_names else NA
      )

    } else {
      # Single-pathway potato
      is_verified <- potato@verified %||% FALSE

      results[[potato_id]] <- list(
        verified = if (is_verified) 1 else 0,
        total = 1,
        unverified = if (!is_verified) potato_id else NA
      )
    }
  }

  # Print summary
  cli::cli_h2("Pathway Verification Status")
  cli::cli_text("")

  for (potato_id in names(results)) {
    res <- results[[potato_id]]

    if (res$verified == res$total) {
      # All verified
      cli::cli_alert_success("{.strong {potato_id}}: {res$verified}/{res$total} verified")
    } else {
      # Some unverified
      cli::cli_alert_warning("{.strong {potato_id}}: {res$verified}/{res$total} verified")
      if (!is.na(res$unverified[1])) {
        cli::cli_text("  {.emph Unverified}: {paste(res$unverified, collapse=', ')}")
      }
    }
  }

  cli::cli_text("")

  # Build tibble for return
  result_df <- dplyr::bind_rows(lapply(names(results), function(id) {
    res <- results[[id]]
    data.frame(
      potato = id,
      verified = res$verified,
      total = res$total,
      unverified = if (!is.na(res$unverified[1])) paste(res$unverified, collapse = ", ") else NA_character_,
      stringsAsFactors = FALSE
    )
  }))

  invisible(result_df)
}
