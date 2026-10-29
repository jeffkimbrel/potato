# Initialize a new potato project

Creates a new project folder with the standard structure: potatoes
directory, config file, RStudio project file, and folders for genomes
and results. A potato project is the on-disk folder structure. Use
create_sack() to load the project into an in-memory PotatoSack object
for analysis.

## Usage

``` r
initialize_potato_project(path, copy_potatoes = TRUE, template_config = NULL)
```

## Arguments

- path:

  Character. Full path where the project will be created.

- copy_potatoes:

  Logical. If TRUE (default), copies bundled potatoes from package to
  project.

- template_config:

  Character. Optional path to existing potato_config.yaml to copy
  database settings from.

## Value

Invisibly returns the path to the new project folder.
