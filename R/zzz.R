.onLoad <- function(libname, pkgname){
  reticulate::use_condaenv("potato")

  reticulate::py_require("jakomics")
  jakomics <<- reticulate::import("jakomics", delay_load = TRUE)

  gator <<- reticulate::import_from_path("metadata", file.path("inst", "python"))

}
