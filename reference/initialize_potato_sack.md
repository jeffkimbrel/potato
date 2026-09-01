# Initialize a new potato sack project

Creates a new project folder with the standard structure: potatoes
directory, config file, RStudio project file, and folders for genomes
and results. A "potato sack" is a collection of potatoes (pathways) for
an annotation project.

## Usage

``` r
initialize_potato_sack(path, copy_potatoes = TRUE)
```

## Arguments

- path:

  Character. Full path where the project will be created.

- copy_potatoes:

  Logical. If TRUE (default), copies example potatoes from package to
  project.

## Value

Invisibly returns the path to the new project folder.
