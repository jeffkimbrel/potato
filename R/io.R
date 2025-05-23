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
#' @returns A formatted potato DB
#' @export

get_potato_db = function(file_path = "potato") {

  # if file_path is not "internal", then use the file_path to create a potato_db
  if (file_path == "potato") {
    potato_db = potato.m$Metadata(system.file("extdata", "potato_db.xlsx", package = "potato"))
    potato_db$summary()
    return(potato_db)
  } else {

    if (!file.exists(file_path)) {
      stop("File does not exist at ", file_path, call.= FALSE)
    }

    potato_db = potato.m$Metadata(file_path)
    potato_db$summary()
    return(potato_db)
  }
}
