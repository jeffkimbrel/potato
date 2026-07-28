#' Initialize annotation results table
#'
#' Creates a tibble for storing genome annotation results across multiple tools.
#' Each row is a genome, with list columns for tool results. Tool columns are
#' inferred from database types in config.
#'
#' @param genomes List of jakomics FILE objects
#' @param config Potato config object (optional). If provided, tool columns are
#'   inferred from database types. If NULL, uses default tools.
#' @param tools Character vector of tool names to create columns for.
#'   Only used if config is NULL. Default: c("kofam", "pfam", "hmm", "blast")
#'
#' @returns tibble with columns:
#'   - genome: genome short name
#'   - file_object: jakomics FILE object
#'   - One column per tool type found in config (empty list initially)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- load_potato_config()
#' genomes <- prepare_genomes(genome_dir = "genomes/")
#' anno_table <- initialize_annotation_table(genomes, config)
#'
#' # Run tools and populate columns
#' anno_table$kofam <- annotate_genomes_kofam(genomes, potatoes, config)
#'
#' # Unnest specific tool results
#' library(tidyr)
#' unnest(anno_table, kofam)
#' }
initialize_annotation_table <- function(genomes, config = NULL,
                                       tools = c("kofam", "pfam", "hmm", "blast")) {

  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required. Install with: install.packages('tibble')", call. = FALSE)
  }

  if (length(genomes) == 0) {
    stop("No genomes provided", call. = FALSE)
  }

  # Infer database names from config
  if (!is.null(config)) {
    if (!is.null(config$databases)) {
      # Use database names, not types
      db_names <- names(config$databases)
      tools <- db_names
      message("Initializing annotation table with databases from config: ", paste(tools, collapse = ", "))
    }
  }

  # Extract genome names
  genome_names <- sapply(genomes, function(g) g$short_name)

  # Build initial tibble
  anno_table <- tibble::tibble(
    genome = genome_names,
    file_object = genomes
  )

  # Add empty list columns for each database
  for (db_name in tools) {
    anno_table[[db_name]] <- vector("list", length(genomes))
  }

  structure(anno_table, class = c("potato_annotation_table", class(anno_table)))
}


#' Print method for annotation table
#' @export
print.potato_annotation_table <- function(x, ...) {
  cat("<Potato Annotation Table>\n")
  cat("  Genomes:", nrow(x), "\n")

  # Count which tools have results
  tool_cols <- setdiff(names(x), c("genome", "file_object"))
  for (tool in tool_cols) {
    n_complete <- sum(sapply(x[[tool]], function(r) !is.null(r) && nrow(r) > 0))
    if (n_complete > 0) {
      cat("  ", tool, ":", n_complete, "of", nrow(x), "genomes annotated\n")
    }
  }

  cat("\n")
  NextMethod()
}


#' Get summary of annotation results
#'
#' Summarizes how many hits each genome has for each tool
#'
#' @param anno_table Annotation table from initialize_annotation_table()
#'
#' @returns tibble with summary statistics
#' @keywords internal
summarize_annotations <- function(anno_table) {
  if (!inherits(anno_table, "potato_annotation_table")) {
    stop("Input must be a potato_annotation_table", call. = FALSE)
  }

  tool_cols <- setdiff(names(anno_table), c("genome", "file_object"))

  summary_list <- lapply(tool_cols, function(tool) {
    counts <- sapply(anno_table[[tool]], function(result) {
      if (is.null(result) || !is.data.frame(result)) return(0)
      nrow(result)
    })

    passed_counts <- sapply(anno_table[[tool]], function(result) {
      if (is.null(result) || !is.data.frame(result) || !"passed" %in% names(result)) return(NA)
      sum(result$passed, na.rm = TRUE)
    })

    tibble::tibble(
      genome = anno_table$genome,
      tool = tool,
      total_hits = counts,
      passed_threshold = passed_counts
    )
  })

  do.call(rbind, summary_list)
}


#' Unnest all tool results into long format
#'
#' Combines all tool results into a single data frame with tool column
#'
#' @param anno_table Annotation table
#' @param tools Character vector of which tools to include. NULL = all tools
#'
#' @returns Long-format data frame with all annotation results
#' @keywords internal
unnest_all_annotations <- function(anno_table, tools = NULL) {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Install with: install.packages('tidyr')", call. = FALSE)
  }

  if (!inherits(anno_table, "potato_annotation_table")) {
    stop("Input must be a potato_annotation_table", call. = FALSE)
  }

  tool_cols <- setdiff(names(anno_table), c("genome", "file_object"))

  if (!is.null(tools)) {
    tool_cols <- intersect(tool_cols, tools)
  }

  # For each tool, unnest and add tool column
  result_list <- lapply(tool_cols, function(tool) {
    # Filter to non-empty results
    has_results <- sapply(anno_table[[tool]], function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0)

    if (!any(has_results)) return(NULL)

    subset_table <- anno_table[has_results, c("genome", tool)]

    # Remove genome column from nested data if present (avoid duplicates)
    subset_table[[tool]] <- lapply(subset_table[[tool]], function(df) {
      if ("genome" %in% names(df)) {
        df[, setdiff(names(df), "genome"), drop = FALSE]
      } else {
        df
      }
    })

    unnested <- tidyr::unnest(subset_table, cols = !!rlang::sym(tool))
    unnested$tool <- tool
    unnested
  })

  # Remove NULLs and combine
  result_list <- result_list[!sapply(result_list, is.null)]

  if (length(result_list) == 0) {
    return(data.frame())
  }

  # Use dplyr::bind_rows to handle different column structures from different tools
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    # Fallback: find common columns
    all_cols <- unique(unlist(lapply(result_list, names)))
    result_list <- lapply(result_list, function(df) {
      missing_cols <- setdiff(all_cols, names(df))
      for (col in missing_cols) {
        df[[col]] <- NA
      }
      df[, all_cols]
    })
    do.call(rbind, result_list)
  } else {
    dplyr::bind_rows(result_list)
  }
}
