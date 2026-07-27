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

    #' @field genomes List of jakomics FILE objects
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


#' Create provenance record
#'
#' Records command, parameters, tool versions, and timestamps for reproducibility.
#'
#' @param stage_name Name of the stage (e.g., "annotate_kofam", "score_pathways")
#' @param command Command string or description
#' @param params List of parameters used
#' @param versions List of tool versions (e.g., list(kofamscan = "1.3.0"))
#' @param database_info List of database info (paths, hashes)
#' @param files Optional list of input/output files
#'
#' @return List with provenance information
#' @importFrom utils packageVersion
#' @keywords internal
create_provenance <- function(stage_name, command, params = list(),
                             versions = list(), database_info = list(),
                             files = list()) {
  prov <- list(
    stage = stage_name,
    command = command,
    params = params,
    versions = versions,
    database_info = database_info,
    timestamp = Sys.time(),
    potato_version = as.character(packageVersion("potato")),
    r_version = R.version.string,
    platform = R.version$platform
  )

  # Add file lists if provided
  if (length(files) > 0) {
    prov$files <- files
  }

  prov
}


#' Get tool version
#'
#' Attempts to get version of external tool by running --version
#'
#' @param tool_name Name of the tool (e.g., "exec_annotation", "blastp")
#' @return Version string or "unknown"
#' @keywords internal
get_tool_version <- function(tool_name) {
  tryCatch({
    # Try common version flags
    for (flag in c("--version", "-version", "-v", "version")) {
      result <- system2(tool_name, flag, stdout = TRUE, stderr = TRUE)
      if (!is.null(attr(result, "status")) && attr(result, "status") != 0) next
      if (length(result) > 0) {
        # Extract version number from first line
        version_line <- result[1]
        return(version_line)
      }
    }
    "unknown"
  }, error = function(e) "unknown")
}


#' Calculate file hash
#'
#' Calculate MD5 hash of a file for change detection
#'
#' @param file_path Path to file
#' @return MD5 hash string or NULL if file doesn't exist
#' @keywords internal
calculate_file_hash <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)

  tryCatch({
    digest::digest(file = file_path, algo = "md5")
  }, error = function(e) NULL)
}


#' Calculate directory hash
#'
#' Calculate combined hash of all files in a directory
#'
#' @param dir_path Path to directory
#' @param pattern Optional file pattern to match
#' @return MD5 hash string or NULL
#' @keywords internal
calculate_dir_hash <- function(dir_path, pattern = NULL) {
  if (!dir.exists(dir_path)) return(NULL)

  tryCatch({
    files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) return(NULL)

    # Hash of all file hashes combined
    file_hashes <- sapply(files, calculate_file_hash)
    digest::digest(paste(file_hashes, collapse = ""), algo = "md5")
  }, error = function(e) NULL)
}
