#' Potato Sack S7 Class
#'
#' An S7 class representing a complete potato annotation project, including
#' configuration, potatoes, genomes, annotation results, and provenance tracking.
#'
#' @export
PotatoSack <- S7::new_class(
  name = "PotatoSack",
  package = "potato",
  properties = list(
    #' @field sack_id Unique identifier for this sack (hash or user-provided)
    sack_id = S7::class_character,

    #' @field sack_root Path to the sack directory
    sack_root = S7::class_character,

    #' @field config Loaded potato configuration
    config = S7::class_any,

    #' @field potatoes List of Potato objects
    potatoes = S7::new_property(S7::class_list, default = list()),

    #' @field genomes List of GenomeFile objects
    genomes = S7::new_property(S7::class_list, default = list()),

    #' @field results Annotation results table (tibble)
    results = S7::new_property(S7::class_any, default = NULL),

    #' @field scores Pathway scoring results (tibble)
    scores = S7::new_property(S7::class_any, default = NULL),

    #' @field completed_stages Vector of completed annotation stages
    completed_stages = S7::new_property(S7::class_character, default = character(0)),

    #' @field provenance Provenance tracking for each stage
    provenance = S7::new_property(S7::class_list, default = list()),

    #' @field metadata Additional metadata
    metadata = S7::new_property(S7::class_list, default = list())
  )
)


#' Print method for PotatoSack
#' @export
S7::method(print, PotatoSack) <- function(x, ...) {
  if (!requireNamespace("cli", quietly = TRUE)) {
    cat("<Potato Sack>\n")
    cat("  ID:", x@sack_id, "\n")
    cat("  Location:", x@sack_root, "\n")
    cat("  Potatoes:", length(x@potatoes), "\n")
    cat("  Genomes:", length(x@genomes), "\n")
    if (length(x@completed_stages) > 0) {
      cat("  Stages:", paste(x@completed_stages, collapse = " -> "), "\n")
    }
  } else {
    cli::cli_h1("Potato Sack")
    cli::cli_dl(c(
      "ID" = x@sack_id,
      "Location" = x@sack_root,
      "Potatoes" = length(x@potatoes),
      "Genomes" = length(x@genomes)
    ))

    if (length(x@completed_stages) > 0) {
      cli::cli_text("Stages: {paste(x@completed_stages, collapse = ' -> ')}")
    }

    if (!is.null(x@config$databases)) {
      cli::cli_h2("Databases")
      for (db_name in names(x@config$databases)) {
        db <- x@config$databases[[db_name]]
        cli::cli_text("  {cli::symbol$bullet} {db_name} ({db$type})")
      }
    }

    if (!is.null(x@results)) {
      cli::cli_h2("Results")
      cli::cli_text("Annotation table: {nrow(x@results)} genome{?s}")
    }

    if (!is.null(x@scores)) {
      cli::cli_h2("Scores")
      cli::cli_text("Pathway scores: {nrow(x@scores)} result{?s}")
    }
  }

  invisible(x)
}


#' Summary method for PotatoSack
#' @export
S7::method(summary, PotatoSack) <- function(object, ...) {
  print(object)

  # Show provenance summary if available
  if (length(object@provenance) > 0 && requireNamespace("cli", quietly = TRUE)) {
    cli::cli_h2("Provenance")
    for (stage in names(object@provenance)) {
      prov <- object@provenance[[stage]]
      cli::cli_text("{cli::symbol$arrow_right} {stage}: {prov$timestamp}")
    }
  }

  invisible(object)
}

