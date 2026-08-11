#' Plot potato network (interactive visNetwork)
#'
#' Creates an interactive visNetwork visualization of a potato pathway network.
#' Supports drag-and-drop node positioning, coordinate export, and genome detection status.
#'
#' @param potato Potato object or path to potato JSON file
#' @param sack PotatoSack object (optional, for genome detection status)
#' @param genome_name Genome name to show detection status (requires sack)
#' @param show_compounds Logical. Show compound nodes in bipartite graph (default: FALSE)
#' @param pathway For multi-pathway networks, show only this pathway (pathway ID). NULL = show all
#' @param layout Layout algorithm. "curated" uses x,y coordinates from potato JSON (default), or specify igraph layout: "fr", "kk", "dh", "circle", "nicely"
#'
#' @return A visNetwork htmlwidget
#' @export
plot_potato_interactive <- function(potato, sack = NULL, genome_name = NULL,
                                     show_compounds = FALSE, pathway = NULL, layout = "curated") {

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg visNetwork} is required for interactive plots",
      "i" = "Install with: install.packages('visNetwork')"
    ))
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg igraph} is required for plotting")
  }

  # Prepare data
  prep <- prepare_potato_for_plotting(potato, sack, genome_name, show_compounds, pathway)

  # Calculate layout (scale up coordinates for visNetwork pixel space)
  node_coords <- calculate_node_layout(prep$potato, prep$g, prep$is_multi_pathway, show_compounds,
                                        layout = layout, scale_for_visnetwork = TRUE)

  # Build visNetwork
  build_visnetwork(prep$potato, prep$g, node_coords, prep$node_status, prep$has_genome)
}
