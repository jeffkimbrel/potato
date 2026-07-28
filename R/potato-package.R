#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats median quantile
## usethis namespace: end
NULL

# Suppress R CMD check notes for NSE in dplyr
utils::globalVariables(c(
  "genome", "potato", "potato_id", "potato_gene", "database", "passed",
  "n", "present", "completeness", "detection_rate", "fill_color"
))

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL
