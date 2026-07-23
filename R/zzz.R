.onLoad <- function(libname, pkgname){
  # Only load jakomics for tool execution
  # Use delay_load so package loads even if conda env not active
  tryCatch({
    reticulate::use_condaenv("potato")
    reticulate::py_require("jakomics")
    jakomics <<- reticulate::import("jakomics", delay_load = TRUE)
  }, error = function(e) {
    packageStartupMessage("Note: jakomics not loaded. Tool execution will be unavailable.")
    packageStartupMessage("To enable annotation: activate 'potato' conda environment")
  })
}
