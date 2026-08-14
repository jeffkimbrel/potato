#' Normalize compound names by sorting multiple compounds
#' @noRd
normalize_compound_name <- function(compound_name) {
  if (is.null(compound_name)) {
    return(compound_name)
  }
  if (is.na(compound_name) || nchar(compound_name) == 0) {
    return(compound_name)
  }

  # Split by " + " (with spaces), sort, rejoin
  parts <- strsplit(compound_name, " \\+ ")[[1]]
  if (length(parts) > 1) {
    parts <- sort(trimws(parts))
    return(paste(parts, collapse = " + "))
  }
  return(compound_name)
}