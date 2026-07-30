#' POTATO plot theme
#'
#' Consistent theme for all POTATO visualizations with transparent backgrounds.
#'
#' @export
potato_theme <- function() {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent")
  )
}
