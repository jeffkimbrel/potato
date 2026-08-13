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
    config_path <- file.path(sack_root, "potato_config.yaml")
  }

  if (!file.exists(config_path)) {
    cli::cli_abort("Config file not found: {.path {config_path}}")
  }
  # Load YAML
  config <- yaml::read_yaml(config_path)

  # Validate required sections
  if (is.null(config$databases)) {
    cli::cli_abort("Config missing {.field databases:} section")
  }

  # Validate and normalize database configs
  config$databases <- validate_database_configs(config$databases)

  # Add config path for reference
  config$config_path <- normalizePath(config_path)

  # returning config, but making it a structure just returns it with a class designation
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
    cli::cli_abort("Config contains duplicate database names: {.val {duplicates}}")
  }

  for (db_name in names(databases)) {
    db <- databases[[db_name]]

    # Check type is specified
    if (is.null(db$type)) {
      cli::cli_abort("Database {.val {db_name}} missing {.field type} field")
    }

    # Validate based on type
    if (db$type == "kofam") {
      # Check profiles_dir and ko_list
      if (!is.null(db$profiles_dir) && !dir.exists(db$profiles_dir)) {
        cli::cli_abort("Database {.val {db_name}}: profiles_dir not found: {.path {db$profiles_dir}}")
      }
      if (!is.null(db$ko_list) && !file.exists(db$ko_list)) {
        cli::cli_abort("Database {.val {db_name}}: ko_list not found: {.path {db$ko_list}}")
      }
      # Set defaults
      if (is.null(db$executable)) databases[[db_name]]$executable <- "exec_annotation"

    } else if (db$type == "blast") {
      # Check path exists
      if (!is.null(db$path) && !file.exists(db$path)) {
        cli::cli_abort("Database {.val {db_name}}: path not found: {.path {db$path}}")
      }
      # Set defaults
      if (is.null(db$executable)) databases[[db_name]]$executable <- "blastp"

    } else if (db$type == "hmm" || db$type == "pfam") {
      # Check path exists
      if (!is.null(db$path) && !file.exists(db$path)) {
        cli::cli_abort("Database {.val {db_name}}: path not found: {.path {db$path}}")
      }
      # Set defaults
      if (is.null(db$executable)) databases[[db_name]]$executable <- "hmmsearch"

    } else {
      cli::cli_abort("Database {.val {db_name}}: unknown type '{db$type}'")
    }
  }

  databases
}
