# Load multiple potatoes from directory

By default, only loads active potatoes (active != false). Set
include_inactive = TRUE to load inactive potatoes.

## Usage

``` r
load_potatoes(dir, tags = NULL, include_inactive = FALSE)
```

## Arguments

- dir:

  Directory containing potato JSON files

- tags:

  Optional character vector of tags to filter by

- include_inactive:

  Logical. If TRUE, loads inactive potatoes (active: false)

## Value

Named list of Potato objects
