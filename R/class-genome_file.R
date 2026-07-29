#' Genome File S7 Class
#'
#' An S7 class representing a genome file with metadata extracted from jakomics FILE objects.
#' This class is serialization-safe (unlike Python objects from jakomics).
#'
#' @export
GenomeFile <- S7::new_class(
  name = "GenomeFile",
  package = "potato",
  properties = list(
    #' @field short_name Short genome name (e.g., "GCA_001234")
    short_name = S7::class_character,

    #' @field file_path Absolute path to genome file
    file_path = S7::class_character,

    #' @field name Full file name with extension
    name = S7::class_character,

    #' @field file_type File extension (e.g., "faa", "fasta")
    file_type = S7::class_character,

    #' @field md5 MD5 hash of file contents (optional, for change tracking)
    md5 = S7::new_property(S7::class_character, default = "")
  )
)


#' Print method for GenomeFile
#' @export
S7::method(print, GenomeFile) <- function(x, ...) {
  cat("<GenomeFile>\n")
  cat("  Name:", x@short_name, "\n")
  cat("  Path:", x@file_path, "\n")
  cat("  Type:", x@file_type, "\n")
  if (nchar(x@md5) > 0) {
    cat("  MD5:", x@md5, "\n")
  }
  invisible(x)
}


#' Convert jakomics FILE object to GenomeFile S7 object
#' @param jakomics_file A jakomics FILE Python object
#' @return GenomeFile S7 object
#' @noRd
jakomics_to_genome_file <- function(jakomics_file) {
  # Extract file extension without the leading dot
  file_type <- sub("^\\.", "", jakomics_file$suffix)

  GenomeFile(
    short_name = jakomics_file$short_name,
    file_path = jakomics_file$file_path,
    name = jakomics_file$name,
    file_type = file_type,
    md5 = ""  # Could compute if needed: jakomics_file$get_md5()
  )
}
