.onLoad <- function(libname, pkgname){
  # Determine which conda environment to use
  # 1. Check environment variable first (allows user override)
  # 2. Default to "potato"
  conda_env <- Sys.getenv("POTATO_CONDA_ENV", unset = "potato")

  # Try to load jakomics for tool execution
  # Suppress errors here - we'll notify in .onAttach if needed
  tryCatch({
    # Force reticulate to use conda environment (required = TRUE)
    # This prevents falling back to system Python
    reticulate::use_condaenv(conda_env, required = TRUE)

    # Verify we're using the right Python
    py_config <- reticulate::py_config()
    if (!grepl("conda|miniconda|miniforge", py_config$python, ignore.case = TRUE)) {
      stop("Not using conda Python. Found: ", py_config$python)
    }

    reticulate::py_require("jakomics")
    # Store in package namespace, not global environment
    assign("jakomics", reticulate::import("jakomics", delay_load = TRUE), envir = parent.env(environment()))
  }, error = function(e) {
    # Silent - will notify in .onAttach
  })
}

.onAttach <- function(libname, pkgname) {
  # Check if jakomics was loaded successfully in package namespace
  jakomics_loaded <- exists("jakomics", envir = asNamespace(pkgname), inherits = FALSE)

  if (!jakomics_loaded) {
    packageStartupMessage("Note: jakomics not loaded. Tool execution will be unavailable.")
    packageStartupMessage("")
    packageStartupMessage("To enable annotation tools:")
    packageStartupMessage("  1. Create conda environment from package:")
    packageStartupMessage("     env_file <- system.file('environment.yaml', package = 'potato')")
    packageStartupMessage("     system(paste('conda env create -f', env_file))")
    packageStartupMessage("     # This creates an environment named 'potato'")
    packageStartupMessage("")
    packageStartupMessage("  2. Restart R session")
    packageStartupMessage("")
    packageStartupMessage("To use a custom environment name:")
    packageStartupMessage("  1. Create with custom name:")
    packageStartupMessage("     system(paste('conda env create -f', env_file, '-n my_custom_name'))")
    packageStartupMessage("  2. Set environment variable before loading potato:")
    packageStartupMessage("     Sys.setenv(POTATO_CONDA_ENV = 'my_custom_name')")
    packageStartupMessage("     # Or add to .Renviron: POTATO_CONDA_ENV=my_custom_name")
  }
}
