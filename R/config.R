#' Load potato configuration from YAML
#'
#' Reads and validates a potato_config.yaml file. Returns a structured list
#' with tool configurations, database paths, and thresholds.
#'
#' @param config_path Path to potato_config.yaml file. If NULL, searches for it
#'   using find_potato_sack() from current directory.
#'
#' @returns List with parsed configuration
#'
#' @details
#' The configuration file defines:
#' - Tool executables and database paths (kofam, pfam, hmmer, blast)
#' - Default thresholds for each tool
#' - Project paths and settings
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- load_potato_config("potato_config.yaml")
#' config <- load_potato_config()  # Auto-finds in current sack
#' }
load_potato_config <- function(config_path = NULL) {

  # Auto-find config if not specified
  if (is.null(config_path)) {
    sack_root <- find_potato_sack()
    if (is.null(sack_root)) {
      stop("Not inside a potato sack. Provide config_path or run from within a sack.", call. = FALSE)
    }
    config_path <- file.path(sack_root, "potato_config.yaml")
  }

  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path, call. = FALSE)
  }

  # Load YAML
  config <- yaml::read_yaml(config_path)

  # Validate required sections
  if (is.null(config$databases)) {
    stop("Config missing 'databases:' section", call. = FALSE)
  }

  # Validate and normalize database configs
  config$databases <- validate_database_configs(config$databases)

  # Add config path for reference
  config$config_path <- normalizePath(config_path)

  structure(config, class = "potato_config")
}


#' Validate database configurations
#'
#' @param databases Databases section from config
#' @keywords internal
validate_database_configs <- function(databases) {

  # Check for duplicate database names
  db_names <- names(databases)
  if (length(db_names) != length(unique(db_names))) {
    duplicates <- db_names[duplicated(db_names)]
    stop("Config contains duplicate database names: ", paste(duplicates, collapse = ", "), call. = FALSE)
  }

  for (db_name in names(databases)) {
    db <- databases[[db_name]]

    # Check type is specified
    if (is.null(db$type)) {
      stop("Database '", db_name, "' missing 'type' field", call. = FALSE)
    }

    # Validate based on type
    if (db$type == "kofam") {
      # Check profiles_dir and ko_list
      if (!is.null(db$profiles_dir) && !dir.exists(db$profiles_dir)) {
        warning("Database '", db_name, "': profiles_dir not found: ", db$profiles_dir)
      }
      if (!is.null(db$ko_list) && !file.exists(db$ko_list)) {
        warning("Database '", db_name, "': ko_list not found: ", db$ko_list)
      }
      # Set defaults
      if (is.null(db$executable)) databases[[db_name]]$executable <- "exec_annotation"

    } else if (db$type == "blast") {
      # Check path exists
      if (!is.null(db$path) && !file.exists(db$path)) {
        warning("Database '", db_name, "': path not found: ", db$path)
      }
      # Set defaults
      if (is.null(db$executable)) databases[[db_name]]$executable <- "blastp"

    } else if (db$type == "hmm" || db$type == "pfam") {
      # Check path exists
      if (!is.null(db$path) && !file.exists(db$path)) {
        warning("Database '", db_name, "': path not found: ", db$path)
      }
      # Set defaults
      if (is.null(db$executable)) databases[[db_name]]$executable <- "hmmsearch"

    } else {
      warning("Database '", db_name, "': unknown type '", db$type, "'")
    }
  }

  databases
}


#' Print method for potato_config
#' @export
print.potato_config <- function(x, ...) {
  cat("<Potato Configuration>\n")
  cat("  Config:", x$config_path, "\n")

  if (!is.null(x$project_name)) {
    cat("  Project:", x$project_name, "\n")
  }

  cat("\nDatabases:\n")

  # Count by type
  db_types <- sapply(x$databases, function(db) db$type)
  type_counts <- table(db_types)

  for (type in names(type_counts)) {
    cat("  ", type, ": ", type_counts[type], " database(s)\n", sep = "")
  }

  cat("\nConfigured databases:\n")
  for (db_name in names(x$databases)) {
    cat("  - ", db_name, " (", x$databases[[db_name]]$type, ")\n", sep = "")
  }

  invisible(x)
}
