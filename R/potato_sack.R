#' Initialize a new potato sack project
#'
#' Creates a new project folder with the standard structure: potatoes directory,
#' config file, RStudio project file, and folders for genomes and results.
#' A "potato sack" is a collection of potatoes (pathways) for an annotation project.
#'
#' @param path Character. Full path where the project will be created.
#' @param copy_potatoes Logical. If TRUE (default), copies example potatoes from package to project.
#'
#' @returns Invisibly returns the path to the new project folder.
#' @export
initialize_potato_sack <- function(path, copy_potatoes = TRUE) {

  project_path <- normalizePath(path, mustWork = FALSE)
  name <- basename(project_path)

  if (dir.exists(project_path)) {
    stop("Directory already exists: ", project_path, "\n",
         "Choose a different name or remove the existing folder.",
         call. = FALSE)
  }

  # --- Folder structure ---
  dirs <- c(
    project_path,
    file.path(project_path, "potatoes"),
    file.path(project_path, "genomes"),
    file.path(project_path, "results"),
    file.path(project_path, "results", "annotations"),
    file.path(project_path, "results", "scores"),
    file.path(project_path, ".claude")
  )

  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # --- Copy example potatoes from package ---
  if (copy_potatoes) {
    package_potatoes <- system.file("potatoes", package = "potato")
    if (nzchar(package_potatoes) && dir.exists(package_potatoes)) {
      potato_files <- list.files(package_potatoes, pattern = "\\.json$", full.names = TRUE)
      if (length(potato_files) > 0) {
        file.copy(potato_files, file.path(project_path, "potatoes"), overwrite = FALSE)
        cat("  Copied", length(potato_files), "example potatoes\n")
      }
    }
  }

  # --- Copy Claude agents from package ---
  package_claude <- system.file(".claude", package = "potato")
  if (nzchar(package_claude) && dir.exists(package_claude)) {
    # Copy entire .claude directory structure
    claude_files <- list.files(package_claude, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    if (length(claude_files) > 0) {
      for (src_file in claude_files) {
        # Get relative path from package .claude dir
        rel_path <- sub(paste0("^", package_claude, "/"), "", src_file)
        dest_file <- file.path(project_path, ".claude", rel_path)

        # Create parent directory if needed
        dest_dir <- dirname(dest_file)
        if (!dir.exists(dest_dir)) {
          dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
        }

        file.copy(src_file, dest_file, overwrite = FALSE)
      }
      cat("  Copied Claude agents to .claude/\n")
    }
  }

  # --- potato_config.yaml ---
  yaml_lines <- c(
    paste0("project_name: ", name),
    paste0("created: ", format(Sys.Date(), "%Y-%m-%d")),
    "",
    "# Directory paths (relative to this file)",
    "paths:",
    "  potatoes: potatoes/",
    "  genomes: genomes/",
    "  results: results/",
    "",
    "# Database definitions",
    "# Each database has a type (kofam, blast, hmm, pfam) and type-specific config",
    "databases:",
    "  # KOfam database",
    "  kofam:",
    "    type: kofam",
    "    profiles_dir: /path/to/kofam/profiles/",
    "    ko_list: /path/to/kofam/ko_list",
    "    executable: exec_annotation",
    "    thresholds:",
    "      score_ratio: 1.0  # 1.0 = use native KOfam thresholds",
    "  ",
    "  # BLAST database",
    "  # blast:",
    "  #   type: blast",
    "  #   path: /path/to/database.faa",
    "  #   executable: blastp",
    "  #   thresholds:",
    "  #     evalue: 1e-7",
    "  #     bitscore: 50",
    "  #     pident: 30",
    "  ",
    "  # HMM database",
    "  # hmm:",
    "  #   type: hmm",
    "  #   path: /path/to/profiles.hmm",
    "  #   executable: hmmsearch",
    "  #   thresholds:",
    "  #     evalue: 1e-10",
    "  ",
    "  # PFAM database",
    "  # pfam:",
    "  #   type: pfam",
    "  #   path: /path/to/Pfam-A.hmm",
    "  #   executable: hmmsearch",
    "  #   thresholds:",
    "  #     evalue: 1e-5",
    "",
    "# Annotation settings",
    "annotation:",
    "  parallel: true",
    "  workers: 4",
    "",
    "# Scoring settings",
    "scoring:",
    "  min_fraction: 0.75",
    "  use_specificity: true"
  )
  writeLines(yaml_lines, file.path(project_path, "potato_config.yaml"))

  # --- .Rproj file ---
  rproj_lines <- c(
    "Version: 1.0",
    "",
    "RestoreWorkspace: Default",
    "SaveWorkspace: Default",
    "AlwaysSaveHistory: Default",
    "",
    "EnableCodeIndexing: Yes",
    "UseSpacesForTab: Yes",
    "NumSpacesForTab: 2",
    "Encoding: UTF-8",
    "",
    "RnwWeave: Sweave",
    "LaTeX: pdfLaTeX"
  )
  writeLines(rproj_lines, file.path(project_path, paste0(name, ".Rproj")))

  # --- .gitignore ---
  gitignore_lines <- c(
    ".DS_Store",
    "*.tmp",
    ".Rproj.user/",
    ".Rhistory",
    ".RData",
    "",
    "# Large result files (uncomment to exclude from git)",
    "# results/",
    "# genomes/",
    "",
    "# Personal database paths (uncomment to avoid committing local paths)",
    "# potato_config.yaml"
  )
  writeLines(gitignore_lines, file.path(project_path, ".gitignore"))

  # --- README.md ---
  readme_lines <- c(
    paste0("# ", name),
    "",
    "A potato sack annotation project.",
    "",
    "## Getting Started",
    "",
    "1. **Edit configuration**: Open `potato_config.yaml` and update the database paths under `tools:`",
    "2. **Add genomes**: Place your `.faa` or `.gbk` files in `genomes/`",
    "3. **Select potatoes**: Review pathways in `potatoes/`, add/remove as needed",
    "4. **Run annotation**: Load the R package and run your annotation workflow",
    "",
    "## Project Structure",
    "",
    "```",
    name,
    "|-- potato_config.yaml      # Configuration (database paths, tool settings)",
    "|-- potatoes/               # Pathway definitions (JSON files)",
    "|-- genomes/                # Input genome files (FAA, GBK)",
    "|-- results/                # Output files",
    "|   |-- annotations/        # Gene-level annotation results",
    "|   +-- scores/             # Pathway-level scoring results",
    "|-- .claude/                # Claude Code agents for building potatoes",
    "+-- README.md               # This file",
    "```",
    "",
    "## Configuration",
    "",
    "Before running annotations, you must configure tool paths in `potato_config.yaml`:",
    "",
    "- **kofam**: Set `profiles_dir` and `ko_list` to your local KofamScan installation",
    "- **blast**: Verify `blastp` is in your PATH (or set full path to executable)",
    "- **hmmer**: Verify `hmmsearch` is in your PATH (or set full path to executable)",
    "",
    "To test if tools are available, activate the `potato` conda environment:",
    "",
    "```bash",
    "conda activate potato",
    "which exec_annotation  # Should show kofamscan path",
    "which blastp           # Should show blast path",
    "which hmmsearch        # Should show hmmer path",
    "```",
    "",
    paste0("Created: ", format(Sys.Date(), "%Y-%m-%d"))
  )
  writeLines(readme_lines, file.path(project_path, "README.md"))

  # --- git init ---
  git_out <- system2("git", c("init", project_path), stdout = TRUE, stderr = TRUE)
  git_ok <- !inherits(git_out, "error") && !any(grepl("^fatal:", git_out))
  if (!git_ok) {
    warning("git init failed -- initialize manually if needed.")
  }

  # --- Summary ---
  cat("\n")
  cat("Potato sack initialized: ", project_path, "\n", sep = "")
  cat("  Config   : ", file.path(project_path, "potato_config.yaml"), "\n", sep = "")
  cat("  RStudio  : ", file.path(project_path, paste0(name, ".Rproj")), "\n", sep = "")
  cat("  Potatoes : ", file.path(project_path, "potatoes/"), "\n", sep = "")
  cat("  Genomes  : ", file.path(project_path, "genomes/"), "\n", sep = "")
  cat("  Results  : ", file.path(project_path, "results/"), "\n", sep = "")
  cat("  Agents   : ", file.path(project_path, ".claude/"), "\n", sep = "")
  cat("\n")
  cat("Next steps:\n")
  cat("  1. Open", paste0(name, ".Rproj"), "in RStudio or load folder in Positron\n")
  cat("  2. Edit potato_config.yaml to set your database paths\n")
  cat("  3. Run: sack <- load_potato_sack()\n")
  cat("  4. Add genomes: sack <- add_genomes(sack, \"path/to/genomes/*.faa\")\n")
  cat("  5. Run annotation: results <- annotate_sack_simple(sack)\n")
  cat("\n")

  invisible(project_path)
}


#' Find the root of a potato sack project
#'
#' Walks up the directory tree looking for a `potato_config.yaml` file.
#' Returns the project root path if found, NULL otherwise.
#'
#' @param path Character. Starting path to search from. Defaults to current working directory.
#'
#' @returns Character path to the project root, or NULL if not inside a potato sack.
#' @export
find_potato_sack <- function(path = NULL) {
  if (is.null(path)) {
    path <- getwd()
  }

  path <- normalizePath(path, mustWork = FALSE)

  while (TRUE) {
    if (file.exists(file.path(path, "potato_config.yaml"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) return(NULL)
    path <- parent
  }
}
