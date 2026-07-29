#' Load a potato sack
#'
#' Constructs a PotatoSack S7 object from project files, or loads a saved RDS.
#'
#' @param path Character. Path to sack directory or RDS file. If NULL, searches upward from current directory.
#'
#' @returns PotatoSack S7 object
#' @export
load_potato_sack <- function(path = NULL) {

  # Determine what we're loading
  if (is.null(path)) {
    # Search upward for potato_config.yaml
    sack_root <- find_potato_sack()
    if (is.null(sack_root)) {
      stop("Not in a potato sack directory.\n",
           "  Run initialize_potato_sack() to create one, or\n",
           "  Provide path to existing sack directory",
           call. = FALSE)
    }
    path <- sack_root
  }

  # Check what type of path we have
  if (file.exists(path) && grepl("\\.rds$", path, ignore.case = TRUE)) {
    # Load saved RDS
    message("Loading saved sack from RDS")
    sack <- readRDS(path)
    if (!S7::S7_inherits(sack, PotatoSack)) {
      stop("RDS file does not contain a PotatoSack object", call. = FALSE)
    }
    message("Loaded: ", length(sack@genomes), " genome(s), ", length(sack@potatoes), " potato(es)")
    return(sack)
  }

  # Otherwise, construct from directory
  if (!dir.exists(path)) {
    stop("Directory not found: ", path, call. = FALSE)
  }

  sack_root <- normalizePath(path)
  config_path <- file.path(sack_root, "potato_config.yaml")

  if (!file.exists(config_path)) {
    stop("No potato_config.yaml found in: ", sack_root, call. = FALSE)
  }

  message("Loading potato sack from: ", basename(sack_root))

  # Load config
  config <- load_potato_config(config_path)

  # Load potatoes
  potatoes_dir <- file.path(sack_root, config$paths$potatoes)
  if (!dir.exists(potatoes_dir)) {
    warning("Potatoes directory not found: ", potatoes_dir)
    potatoes <- list()
  } else {
    potatoes <- load_potatoes(potatoes_dir)
    message("  Loaded ", length(potatoes), " potato(es)")
  }

  # Generate sack ID
  sack_id <- basename(sack_root)

  # Create PotatoSack object (genomes are empty - use add_genomes())
  sack <- PotatoSack(
    sack_id = sack_id,
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

  message("  Sack loaded. Use add_genomes() to register genome files.")

  sack
}
