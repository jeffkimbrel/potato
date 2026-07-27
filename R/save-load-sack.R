#' Save a PotatoSack object to disk
#'
#' Saves a PotatoSack object as an RDS file with all metadata, results, and
#' provenance information preserved.
#'
#' @param sack PotatoSack object to save
#' @param path Path to save the RDS file (default: "<sack_id>_potato_sack.rds")
#' @param compress Compression method (default: "gzip")
#'
#' @return Invisibly returns the file path
#' @export
save_potato_sack <- function(sack, path = NULL, compress = "gzip") {

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("Input must be a PotatoSack object", call. = FALSE)
  }

  # Generate default path if not provided
  if (is.null(path)) {
    path <- paste0(sack@sack_id, "_potato_sack.rds")
  }

  # Update save timestamp in metadata
  sack@metadata$last_saved <- Sys.time()

  # Save as RDS
  saveRDS(sack, path, compress = compress)

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Saved PotatoSack to: ", path)
  } else {
    cli::cli_alert_success("Saved PotatoSack to: {.path {path}}")
  }

  invisible(path)
}


#' Load a PotatoSack object from disk
#'
#' Loads a previously saved PotatoSack object from an RDS file.
#'
#' @param path Path to the RDS file
#' @param validate Logical. Validate the loaded object (default TRUE)
#'
#' @return PotatoSack object
#' @export
load_saved_sack <- function(path, validate = TRUE) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  # Load RDS
  sack <- readRDS(path)

  # Validate
  if (validate) {
    if (!S7::S7_inherits(sack, PotatoSack)) {
      stop("Loaded object is not a PotatoSack", call. = FALSE)
    }
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Loaded PotatoSack from: ", path)
  } else {
    cli::cli_alert_success("Loaded PotatoSack: {.val {sack@sack_id}}")
  }

  sack
}


#' Print metadata about a saved PotatoSack
#'
#' Shows information about a saved sack without fully loading it.
#'
#' @param path Path to the RDS file
#'
#' @return Invisibly returns the metadata
#' @export
inspect_saved_sack <- function(path) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  sack <- readRDS(path)

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("File is not a PotatoSack object", call. = FALSE)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    cat("<Saved PotatoSack>\n")
    cat("  File:", path, "\n")
    cat("  ID:", sack@sack_id, "\n")
    cat("  Created:", format(sack@metadata$created), "\n")
    if (!is.null(sack@metadata$last_saved)) {
      cat("  Last saved:", format(sack@metadata$last_saved), "\n")
    }
    cat("  Potatoes:", length(sack@potatoes), "\n")
    cat("  Genomes:", length(sack@genomes), "\n")
    if (length(sack@completed_stages) > 0) {
      cat("  Stages:", paste(sack@completed_stages, collapse = ", "), "\n")
    }
  } else {
    cli::cli_h2("Saved PotatoSack: {.path {basename(path)}}")
    cli::cli_dl(c(
      "ID" = sack@sack_id,
      "Created" = format(sack@metadata$created),
      "Last saved" = if (!is.null(sack@metadata$last_saved)) format(sack@metadata$last_saved) else "N/A",
      "Potatoes" = length(sack@potatoes),
      "Genomes" = length(sack@genomes),
      "Stages" = if (length(sack@completed_stages) > 0) paste(sack@completed_stages, collapse = ", ") else "none"
    ))

    if (!is.null(sack@results)) {
      cli::cli_alert_info("Has annotation results: {nrow(sack@results)} genome{?s}")
    }

    if (!is.null(sack@scores)) {
      cli::cli_alert_info("Has pathway scores: {nrow(sack@scores)} result{?s}")
    }
  }

  invisible(sack@metadata)
}


#' Combine multiple PotatoSacks
#'
#' Merges results from multiple annotation runs. Useful for combining batches
#' of genomes annotated at different times or with different database versions.
#'
#' @param ... PotatoSack objects or paths to saved sacks
#' @param sack_ids Optional character vector of IDs to assign to each sack
#'
#' @return New PotatoSack object with combined results
#' @export
combine_sacks <- function(..., sack_ids = NULL) {

  sacks <- list(...)

  # Load any paths
  sacks <- lapply(sacks, function(s) {
    if (is.character(s)) {
      load_saved_sack(s)
    } else {
      s
    }
  })

  # Validate all are PotatoSacks
  for (i in seq_along(sacks)) {
    if (!S7::S7_inherits(sacks[[i]], PotatoSack)) {
      stop("Argument ", i, " is not a PotatoSack object", call. = FALSE)
    }
  }

  if (length(sacks) < 2) {
    stop("Must provide at least 2 sacks to combine", call. = FALSE)
  }

  # Generate or validate sack_ids
  if (is.null(sack_ids)) {
    sack_ids <- sapply(sacks, function(s) s@sack_id)
  } else if (length(sack_ids) != length(sacks)) {
    stop("sack_ids must have same length as number of sacks", call. = FALSE)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Combining ", length(sacks), " sacks...")
  } else {
    cli::cli_h2("Combining {length(sacks)} sack{?s}")
  }

  # Combine results tables
  combined_results <- NULL
  if (all(sapply(sacks, function(s) !is.null(s@results)))) {
    results_list <- list()

    for (i in seq_along(sacks)) {
      result <- sacks[[i]]@results
      # Add sack_id column
      result$sack_id <- sack_ids[i]
      results_list[[i]] <- result
    }

    combined_results <- do.call(rbind, results_list)

    if (!requireNamespace("cli", quietly = TRUE)) {
      message("  Combined results: ", nrow(combined_results), " genomes")
    } else {
      cli::cli_alert_success("Combined results: {nrow(combined_results)} genome{?s}")
    }
  }

  # Combine scores if available
  combined_scores <- NULL
  if (all(sapply(sacks, function(s) !is.null(s@scores)))) {
    scores_list <- list()

    for (i in seq_along(sacks)) {
      score <- sacks[[i]]@scores
      score$sack_id <- sack_ids[i]
      scores_list[[i]] <- score
    }

    combined_scores <- do.call(rbind, scores_list)

    if (!requireNamespace("cli", quietly = TRUE)) {
      message("  Combined scores: ", nrow(combined_scores), " results")
    } else {
      cli::cli_alert_success("Combined scores: {nrow(combined_scores)} result{?s}")
    }
  }

  # Collect all unique potatoes
  all_potatoes <- list()
  for (sack in sacks) {
    for (potato in sack@potatoes) {
      if (!potato@id %in% names(all_potatoes)) {
        all_potatoes[[potato@id]] <- potato
      }
    }
  }

  # Collect all genomes
  all_genomes <- list()
  for (sack in sacks) {
    all_genomes <- c(all_genomes, sack@genomes)
  }

  # Combine metadata
  combined_metadata <- list(
    combined_from = sack_ids,
    combined_at = Sys.time(),
    source_sacks = lapply(sacks, function(s) {
      list(
        sack_id = s@sack_id,
        created = s@metadata$created,
        databases = s@metadata$databases,
        n_genomes = length(s@genomes),
        n_potatoes = length(s@potatoes)
      )
    })
  )

  # Combine provenance
  combined_provenance <- list()
  for (i in seq_along(sacks)) {
    sack_prov <- sacks[[i]]@provenance
    for (stage in names(sack_prov)) {
      stage_key <- paste0(sack_ids[i], "_", stage)
      combined_provenance[[stage_key]] <- sack_prov[[stage]]
    }
  }

  # Create combined sack
  combined_id <- digest::digest(paste(sack_ids, collapse = "_"), algo = "md5")
  combined_id <- paste0("combined_", substr(combined_id, 1, 8))

  # Use first sack's config as base (warn if they differ)
  base_config <- sacks[[1]]@config

  # Check if configs differ significantly
  config_hashes <- sapply(sacks, function(s) {
    digest::digest(s@config$databases, algo = "md5")
  })

  if (length(unique(config_hashes)) > 1) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      warning("Sacks have different database configurations. Using config from first sack.")
    } else {
      cli::cli_alert_warning("Sacks have different database configurations. Using config from first sack.")
    }
  }

  combined_sack <- PotatoSack(
    sack_id = combined_id,
    sack_root = "(combined)",
    config = base_config,
    potatoes = all_potatoes,
    genomes = all_genomes,
    results = combined_results,
    scores = combined_scores,
    completed_stages = unique(unlist(lapply(sacks, function(s) s@completed_stages))),
    provenance = combined_provenance,
    metadata = combined_metadata
  )

  combined_sack
}
