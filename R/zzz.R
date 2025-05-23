.onLoad <- function(libname, pkgname){
  reticulate::use_condaenv("potato")

  reticulate::py_require("jakomics")
  jakomics <<- reticulate::import("jakomics", delay_load = TRUE)

  potato.m <<- reticulate::import_from_path("metadata", file.path("inst", "python"))
  potato.p <<- reticulate::import_from_path("pathway", file.path("inst", "python"))

}
