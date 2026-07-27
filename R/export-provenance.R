#' Export provenance to executable R script
#'
#' Creates an executable R script that reproduces the potato workflow using
#' the exact commands and parameters from provenance records.
#'
#' @param sack PotatoSack object
#' @param output_file Path to write the R script (default: "<sack_id>_provenance.R")
#' @param include_setup Logical. Include library loading and setup code (default TRUE)
#'
#' @return Invisibly returns the file path
#' @export
export_provenance_script <- function(sack, output_file = NULL, include_setup = TRUE) {

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("Input must be a PotatoSack object", call. = FALSE)
  }

  if (is.null(output_file)) {
    output_file <- paste0(sack@sack_id, "_provenance.R")
  }

  # Start building script
  script <- c(
    "# Potato Provenance Script",
    paste0("# Generated: ", Sys.time()),
    paste0("# Sack ID: ", sack@sack_id),
    "",
    "# =============================================================================="
  )

  if (include_setup) {
    script <- c(script,
      "# SETUP",
      "# ==============================================================================",
      "",
      "library(potato)",
      ""
    )

    # Add conda environment info if available
    if (!is.null(sack@metadata$r_session)) {
      script <- c(script,
        paste0("# R version: ", sack@metadata$r_session$r_version),
        paste0("# Platform: ", sack@metadata$r_session$platform),
        ""
      )
    }
  }

  script <- c(script,
    "# ==============================================================================",
    "# USER CONFIGURATION - Edit paths for your system",
    "# ==============================================================================",
    "",
    paste0("sack_root <- \"", sack@sack_root, "\""),
    ""
  )

  # Add database paths
  if (!is.null(sack@config$databases)) {
    script <- c(script,
      "# Database paths",
      ""
    )

    for (db_name in names(sack@config$databases)) {
      db <- sack@config$databases[[db_name]]

      if (db$type == "kofam") {
        script <- c(script,
          paste0("# ", db_name, " (kofam)"),
          paste0("profiles_dir_", db_name, " <- \"", db$profiles_dir, "\""),
          paste0("ko_list_", db_name, " <- \"", db$ko_list, "\""),
          ""
        )
      } else if (db$type == "blast") {
        script <- c(script,
          paste0("# ", db_name, " (blast)"),
          paste0("blast_db_", db_name, " <- \"", db$path, "\""),
          ""
        )
      } else if (db$type == "hmm" || db$type == "pfam") {
        script <- c(script,
          paste0("# ", db_name, " (", db$type, ")"),
          paste0("hmm_db_", db_name, " <- \"", db$path, "\""),
          ""
        )
      }
    }
  }

  # Add genome paths
  if (length(sack@genomes) > 0) {
    script <- c(script,
      "# Genome files",
      ""
    )

    for (i in seq_along(sack@genomes)) {
      genome <- sack@genomes[[i]]
      script <- c(script,
        paste0("genome_", i, " <- \"", genome$file_path, "\"")
      )
    }

    script <- c(script, "")
  }

  # Add potato paths
  if (length(sack@potatoes) > 0) {
    script <- c(script,
      "# Potato files",
      ""
    )

    for (potato in sack@potatoes) {
      script <- c(script,
        paste0("potato_", potato@id, " <- \"", potato@json_path, "\"")
      )
    }

    script <- c(script, "")
  }

  script <- c(script,
    "# ==============================================================================",
    "# WORKFLOW",
    "# ==============================================================================",
    ""
  )

  # Add workflow steps
  if ("annotation" %in% sack@completed_stages) {
    script <- c(script,
      "# Load sack",
      paste0("sack <- load_potato_sack(\"", sack@sack_root, "\")"),
      "",
      "# Run annotation",
      "sack <- annotate_sack(sack)",
      ""
    )

    # Add tool-specific details from provenance if available
    if (!is.null(sack@provenance$kofam)) {
      prov <- sack@provenance$kofam
      script <- c(script,
        "# KOfam annotation details:",
        paste0("#   Tool version: ", prov$versions$kofamscan %||% "unknown"),
        paste0("#   Timestamp: ", prov$timestamp),
        ""
      )
    }

    if (!is.null(sack@provenance$blast)) {
      prov <- sack@provenance$blast
      script <- c(script,
        "# BLAST annotation details:",
        paste0("#   Tool version: ", prov$versions$blastp %||% "unknown"),
        paste0("#   Timestamp: ", prov$timestamp),
        ""
      )
    }
  }

  if ("scoring" %in% sack@completed_stages) {
    script <- c(script,
      "# Score pathways",
      "sack <- score_pathways(sack)",
      ""
    )
  }

  # Add save step
  script <- c(script,
    "# Save results",
    paste0("save_potato_sack(sack, \"", basename(output_file), ".rds\")"),
    ""
  )

  script <- c(script,
    "# ==============================================================================",
    "# DONE",
    "# ==============================================================================",
    "",
    "cat(\"Provenance script completed\\n\")",
    "print(sack)"
  )

  # Write script
  writeLines(script, output_file)

  if (!requireNamespace("cli", quietly = TRUE)) {
    message("Provenance script written to: ", output_file)
  } else {
    cli::cli_alert_success("Provenance script written to: {.path {output_file}}")
  }

  invisible(output_file)
}


#' View provenance for a specific stage
#'
#' Display detailed provenance information for a specific annotation stage.
#'
#' @param sack PotatoSack object
#' @param stage Stage name (e.g., "annotation", "kofam", "blast", "scoring")
#'
#' @return Invisibly returns the provenance record
#' @export
view_provenance <- function(sack, stage) {

  if (!S7::S7_inherits(sack, PotatoSack)) {
    stop("Input must be a PotatoSack object", call. = FALSE)
  }

  if (!stage %in% names(sack@provenance)) {
    stop("Stage '", stage, "' not found in provenance. Available: ",
         paste(names(sack@provenance), collapse = ", "), call. = FALSE)
  }

  prov <- sack@provenance[[stage]]

  if (!requireNamespace("cli", quietly = TRUE)) {
    cat("=== Provenance: ", stage, " ===\n")
    cat("Timestamp:", format(prov$timestamp), "\n")
    cat("Potato version:", prov$potato_version, "\n\n")

    if (length(prov$versions) > 0) {
      cat("Tool versions:\n")
      for (tool in names(prov$versions)) {
        cat("  ", tool, ": ", prov$versions[[tool]], "\n", sep = "")
      }
      cat("\n")
    }

    if (length(prov$params) > 0) {
      cat("Parameters:\n")
      for (param in names(prov$params)) {
        cat("  ", param, ": ", prov$params[[param]], "\n", sep = "")
      }
      cat("\n")
    }

    if (length(prov$database_info) > 0) {
      cat("Databases:\n")
      for (db in names(prov$database_info)) {
        cat("  ", db, "\n", sep = "")
      }
    }
  } else {
    cli::cli_h2("Provenance: {stage}")
    cli::cli_dl(c(
      "Timestamp" = format(prov$timestamp),
      "Potato version" = prov$potato_version,
      "R version" = prov$r_version
    ))

    if (length(prov$versions) > 0) {
      cli::cli_h3("Tool versions")
      for (tool in names(prov$versions)) {
        cli::cli_text("{cli::symbol$bullet} {tool}: {prov$versions[[tool]]}")
      }
    }

    if (length(prov$params) > 0) {
      cli::cli_h3("Parameters")
      for (param in names(prov$params)) {
        cli::cli_text("{cli::symbol$bullet} {param}: {prov$params[[param]]}")
      }
    }

    if (length(prov$database_info) > 0) {
      cli::cli_h3("Databases")
      for (db in names(prov$database_info)) {
        cli::cli_text("{cli::symbol$bullet} {db}")
      }
    }
  }

  invisible(prov)
}
