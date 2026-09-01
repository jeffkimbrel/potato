#' Get simple pathway verification status table
#'
#' Returns a clean tibble with one row per pathway showing verification status.
#'
#' @param input Either a PotatoSack object or a directory path containing potato JSON files
#'
#' @return Tibble with columns: potato, pathway, verified
#' @export
get_pathway_verification <- function(input) {

  results <- list()

  # Handle PotatoSack
  if (inherits(input, "potato::PotatoSack")) {
    potatoes <- S7::prop(input, "potatoes")

    for (potato in potatoes) {
      for (pathway_id in names(potato@pathways)) {
        pathway <- potato@pathways[[pathway_id]]

        results[[length(results) + 1]] <- list(
          potato = potato@id,
          pathway = pathway_id,
          verified = !is.null(pathway$verified) && pathway$verified == TRUE
        )
      }
    }

  } else if (is.character(input) && length(input) == 1) {
    # Handle directory path
    if (!dir.exists(input)) {
      cli::cli_abort("Directory not found: {input}")
    }

    files <- list.files(input, pattern = "\\.json$", full.names = TRUE)

    if (length(files) == 0) {
      cli::cli_abort("No JSON files found in {input}")
    }

    for (file in files) {
      potato <- load_potato(file)

      for (pathway_id in names(potato@pathways)) {
        pathway <- potato@pathways[[pathway_id]]

        results[[length(results) + 1]] <- list(
          potato = potato@id,
          pathway = pathway_id,
          verified = !is.null(pathway$verified) && pathway$verified == TRUE
        )
      }
    }

  } else {
    cli::cli_abort("Input must be a PotatoSack object or a directory path")
  }

  dplyr::bind_rows(results)
}
