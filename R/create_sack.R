#' Create a PotatoSack from project directory
#'
#' Constructs a PotatoSack S7 object from a project directory containing
#' potato_config.yaml and potatoes. To save/load sacks, use standard R
#' functions: saveRDS(sack, "sack.rds") and readRDS("sack.rds").
#'
#' @param path Character. Path to sack directory. If NULL, searches upward from current directory.
#'
#' @returns PotatoSack S7 object
#' @export

create_sack <- function(path = NULL) {

  # Auto-find sack directory if not specified
  if (is.null(path)) {
    path <- find_potato_sack()
  }

  # Validate directory exists
  if (!dir.exists(path)) {
    cli::cli_abort(c(
      "No potato project found at: {.path {path}}",
      "i" = "Directory does not exist",
      "i" = "Create a new project with: {.code initialize_potato_project('{path}')}"
    ))
  }

  sack_root <- normalizePath(path)
  config_path <- file.path(sack_root, "potato_config.yaml")

  if (!file.exists(config_path)) {
    cli::cli_abort(c(
      "No potato project found at: {.path {sack_root}}",
      "x" = "Missing required file: {.file potato_config.yaml}",
      "i" = "Create a new project with: {.code initialize_potato_project('{sack_root}')}"
    ))
  }

  cli::cli_inform("Creating PotatoSack from: {.path {basename(sack_root)}}")

  # Load config
  config <- load_potato_config(config_path)

  # Load potatoes
  potatoes_dir <- file.path(sack_root, config$paths$potatoes)
  if (!dir.exists(potatoes_dir)) {
    cli::cli_abort(c(
      "Potatoes directory not found: {.path {potatoes_dir}}",
      "i" = "The project structure is incomplete",
      "i" = "Reinitialize with: {.code initialize_potato_project('{sack_root}')}"
    ))
  } else {
    potatoes <- load_potatoes(potatoes_dir) # TODO abilty to utilize tags
    cli::cli_inform("Loaded {length(potatoes)} potato(es)")
  }

  # Create PotatoSack object (genomes are empty - use add_genomes())
  sack <- PotatoSack(
    sack_id = basename(sack_root),
    sack_root = sack_root,
    config = config,
    potatoes = potatoes,
    genomes = list(),
    results = NULL,
    scores = NULL,
    completed_stages = character(0),
    provenance = list(),
    metadata = list(
      created = Sys.time(),
      potato_version = as.character(utils::packageVersion("potato"))
    )
  )

  cli::cli_inform("  Sack created. Use add_genomes() to register genome files.")

  sack
}
