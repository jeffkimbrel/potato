#' Find conda executable (internal)
#' @noRd
find_conda <- function() {
  # 1. Check if conda is already in PATH
  conda_path <- Sys.which("conda")
  if (conda_path != "") {
    return(conda_path)
  }

  # 2. Check CONDA_EXE environment variable (set by conda init)
  conda_exe <- Sys.getenv("CONDA_EXE", unset = "")
  if (conda_exe != "" && file.exists(conda_exe)) {
    return(conda_exe)
  }

  # 3. Search common conda installation locations
  home <- Sys.getenv("HOME")
  common_paths <- c(
    file.path(home, "miniforge3", "bin", "conda"),
    file.path(home, "mambaforge", "bin", "conda"),
    file.path(home, "anaconda3", "bin", "conda"),
    file.path(home, "miniconda3", "bin", "conda"),
    "/opt/anaconda3/bin/conda",
    "/opt/miniconda3/bin/conda",
    "/usr/local/anaconda3/bin/conda",
    "/usr/local/miniconda3/bin/conda"
  )

  for (path in common_paths) {
    if (file.exists(path)) {
      return(path)
    }
  }

  # 4. Not found
  return("")
}
