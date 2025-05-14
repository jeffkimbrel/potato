#' Title
#'
#' @returns
#' @export
#'
#' @examples

get_files = function(file_path) {
  file_list = jakomics$utilities$get_files(
    "",
    file_path,
    c("faa", "gb", "gbk", "gbff"))

  names(file_list) = lapply(file_list, function(x) x$short_name) |> unlist()


  tibble::enframe(file_list, value = "file", name = "short_name")
}
