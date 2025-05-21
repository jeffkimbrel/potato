#' Title
#'
#' @returns A tibble of jakomics file objects
#' @export

get_files = function(file_path) {
  file_list = jakomics$utilities$get_files(
    "",
    file_path,
    c("faa", "gb", "gbk", "gbff"))

  names(file_list) = lapply(file_list, function(x) x$short_name) |> unlist()


  tibble::enframe(file_list, value = "file", name = "short_name")
}




#' Title
#'
#' @param file_path
#'
#' @returns A formatted gator DB
#' @export

get_gator_db = function(file_path = "potato") {

  # if file_path is not "internal", then use the file_path to create a gator_db
  if (file_path == "potato") {
    gator_db = gator$Metadata(system.file("extdata", "gator_db.xlsx", package = "potato"))
    gator_db$summary()
    return(gator_db)
  }
}
