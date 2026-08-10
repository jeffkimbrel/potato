#' Check if potato or pathway is verified before editing
#'
#' Checks verification status and blocks edits to verified pathways unless
#' force=TRUE. This prevents accidental modification of validated pathways.
#'
#' @param potato_path Path to potato JSON file
#' @param pathway_id Optional pathway ID within multi-pathway network
#' @param force Logical. If TRUE, bypass verification check (default: FALSE)
#'
#' @return Logical. TRUE if safe to edit, FALSE (with error) if verified
#' @export
check_verified_status <- function(potato_path, pathway_id = NULL, force = FALSE) {

  if (force) {
    cli::cli_alert_warning("Bypassing verification check (force = TRUE)")
    return(TRUE)
  }

  if (!file.exists(potato_path)) {
    cli::cli_abort("Potato file not found: {potato_path}")
  }

  data <- jsonlite::read_json(potato_path, simplifyVector = FALSE)

  # Check if multi-pathway network
  is_network <- !is.null(data$pathways) && is.list(data$pathways)

  if (is_network) {
    # Multi-pathway: check specific pathway if provided
    if (!is.null(pathway_id)) {
      if (!pathway_id %in% names(data$pathways)) {
        cli::cli_abort("Pathway '{pathway_id}' not found in {basename(potato_path)}")
      }

      pathway <- data$pathways[[pathway_id]]
      if (!is.null(pathway$verified) && pathway$verified == TRUE) {
        cli::cli_abort(c(
          "x" = "Cannot edit verified pathway: {pathway_id}",
          "i" = "This pathway has been manually validated",
          "i" = "To proceed: Set {.code verified = FALSE} in the JSON first, or use {.code force = TRUE}"
        ))
      }
    } else {
      # Check if ANY pathway is verified
      verified_pathways <- character(0)
      for (pw_id in names(data$pathways)) {
        pw <- data$pathways[[pw_id]]
        if (!is.null(pw$verified) && pw$verified == TRUE) {
          verified_pathways <- c(verified_pathways, pw_id)
        }
      }

      if (length(verified_pathways) > 0) {
        cli::cli_alert_warning(c(
          "!" = "This potato contains verified pathways: {.val {verified_pathways}}",
          "i" = "Editing the potato may affect these pathways",
          "i" = "Consider specifying {.arg pathway_id} to check individual pathways"
        ))
      }
    }
  } else {
    # Single-pathway (legacy or inactive): check top-level verified
    if (!is.null(data$verified) && data$verified == TRUE) {
      cli::cli_abort(c(
        "x" = "Cannot edit verified potato: {data$id}",
        "i" = "This potato has been manually validated",
        "i" = "To proceed: Set {.code verified = FALSE} in the JSON first, or use {.code force = TRUE}"
      ))
    }
  }

  return(TRUE)
}


#' Get verification status of potato or pathway
#'
#' Returns verification status without blocking. Useful for reporting.
#'
#' @param potato_path Path to potato JSON file
#' @param pathway_id Optional pathway ID within multi-pathway network
#'
#' @return List with verification details
#' @export
get_verification_status <- function(potato_path, pathway_id = NULL) {

  if (!file.exists(potato_path)) {
    cli::cli_abort("Potato file not found: {potato_path}")
  }

  data <- jsonlite::read_json(potato_path, simplifyVector = FALSE)

  is_network <- !is.null(data$pathways) && is.list(data$pathways)

  if (is_network) {
    if (!is.null(pathway_id)) {
      # Single pathway status
      if (!pathway_id %in% names(data$pathways)) {
        cli::cli_abort("Pathway '{pathway_id}' not found in {basename(potato_path)}")
      }

      pathway <- data$pathways[[pathway_id]]
      verified <- !is.null(pathway$verified) && pathway$verified == TRUE

      return(list(
        potato_id = data$id,
        is_network = TRUE,
        pathway_id = pathway_id,
        verified = verified,
        message = if (verified) "Pathway is verified" else "Pathway is not verified"
      ))
    } else {
      # All pathways status
      statuses <- list()
      for (pw_id in names(data$pathways)) {
        pw <- data$pathways[[pw_id]]
        verified <- !is.null(pw$verified) && pw$verified == TRUE
        statuses[[pw_id]] <- verified
      }

      return(list(
        potato_id = data$id,
        is_network = TRUE,
        pathways = statuses,
        verified_count = sum(unlist(statuses)),
        total_count = length(statuses)
      ))
    }
  } else {
    # Single pathway
    verified <- !is.null(data$verified) && data$verified == TRUE

    return(list(
      potato_id = data$id,
      is_network = FALSE,
      verified = verified,
      message = if (verified) "Potato is verified" else "Potato is not verified"
    ))
  }
}
