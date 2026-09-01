#' Get verification status of pathways
#'
#' Shows which pathways in potatoes are verified or need verification.
#' Can check a single potato or all potatoes in a sack. Optionally blocks
#' editing of verified pathways (guard mode).
#'
#' @param potato_or_sack Potato object, potato path, PotatoSack object, or list of potatoes
#' @param pathway_id Optional pathway ID to check within a potato
#' @param abort_if_verified Logical. If TRUE, aborts with error if pathway is verified (guard mode). Default FALSE.
#' @param force Logical. If TRUE, bypasses verification guard (only used when abort_if_verified=TRUE). Default FALSE.
#'
#' @return Tibble with verification status (or TRUE in guard mode)
#' @export
get_verification_status <- function(potato_or_sack,
                                    pathway_id = NULL,
                                    abort_if_verified = FALSE,
                                    force = FALSE) {

  # Guard mode: bypass check if force=TRUE
  if (abort_if_verified && force) {
    cli::cli_alert_warning("Bypassing verification check (force = TRUE)")
    return(invisible(TRUE))
  }

  # Determine if sack, list of potatoes, or single potato
  if (inherits(potato_or_sack, "potato::PotatoSack")) {
    # PotatoSack object
    potatoes <- S7::prop(potato_or_sack, "potatoes")
  } else if (is.list(potato_or_sack) && length(potato_or_sack) > 0 &&
             !is.null(names(potato_or_sack)) &&
             S7::S7_inherits(potato_or_sack[[1]], PotatoV2)) {
    # Named list of potato objects
    potatoes <- potato_or_sack
  } else {
    # Single potato (object or path)
    if (is.character(potato_or_sack)) {
      if (!file.exists(potato_or_sack)) {
        cli::cli_abort("Potato file not found: {potato_or_sack}")
      }
      potato <- load_potato(potato_or_sack)
    } else {
      potato <- potato_or_sack
    }
    potatoes <- list(potato)
    names(potatoes) <- potato@id
  }

  # If specific pathway requested and only one potato, filter to that pathway
  single_pathway_check <- !is.null(pathway_id) && length(potatoes) == 1

  # Collect results
  results <- list()

  for (potato_id in names(potatoes)) {
    potato <- potatoes[[potato_id]]

    # V2 potatoes always have pathways field
    if (!is.null(potato@pathways) && is.list(potato@pathways)) {

      # If specific pathway requested for this potato
      if (single_pathway_check) {
        if (!pathway_id %in% names(potato@pathways)) {
          cli::cli_abort("Pathway '{pathway_id}' not found in {potato_id}")
        }

        pathway <- potato@pathways[[pathway_id]]
        is_verified <- !is.null(pathway$verified) && pathway$verified == TRUE

        # Guard mode: abort if verified
        if (abort_if_verified && is_verified) {
          cli::cli_abort(c(
            "x" = "Cannot edit verified pathway: {pathway_id}",
            "i" = "This pathway has been manually validated",
            "i" = "To proceed: Set {.code verified = FALSE} in the JSON first, or use {.code force = TRUE}"
          ))
        }

        results[[potato_id]] <- list(
          verified = if (is_verified) 1 else 0,
          total = 1,
          unverified = if (!is_verified) pathway_id else NA
        )
      } else {
        # Check all pathways
        pathway_status <- sapply(names(potato@pathways), function(pw_id) {
          potato@pathways[[pw_id]]$verified %||% FALSE
        })

        n_verified <- sum(pathway_status)
        n_total <- length(pathway_status)
        verified_names <- names(pathway_status)[pathway_status]
        unverified_names <- names(pathway_status)[!pathway_status]

        # Guard mode: warn if ANY pathway is verified (unless specific pathway was requested)
        if (abort_if_verified && n_verified > 0 && is.null(pathway_id)) {
          cli::cli_alert_warning(c(
            "!" = "This potato contains verified pathways: {.val {verified_names}}",
            "i" = "Editing the potato may affect these pathways",
            "i" = "Consider specifying {.arg pathway_id} to check individual pathways"
          ))
        }

        results[[potato_id]] <- list(
          verified = n_verified,
          total = n_total,
          verified_names = if (length(verified_names) > 0) verified_names else NA,
          unverified = if (length(unverified_names) > 0) unverified_names else NA
        )
      }

    } else {
      # Malformed potato - no pathways
      cli::cli_alert_warning("Potato {potato_id} has no pathways field")
      results[[potato_id]] <- list(
        verified = 0,
        total = 0,
        unverified = NA
      )
    }
  }

  # Build tibble for return
  result_tbl <- tibble::tibble(
    potato = names(results),
    verified = sapply(results, function(x) x$verified),
    total = sapply(results, function(x) x$total),
    verified_pathways = sapply(results, function(x) {
      if (!is.na(x$verified_names[1])) paste(x$verified_names, collapse = ", ") else NA_character_
    }),
    unverified = sapply(results, function(x) {
      if (!is.na(x$unverified[1])) paste(x$unverified, collapse = ", ") else NA_character_
    })
  )

  # In guard mode, return TRUE if we got here (no abort)
  if (abort_if_verified) {
    return(invisible(TRUE))
  }

  result_tbl
}
