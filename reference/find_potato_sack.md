# Find the root of a potato sack project

Walks up the directory tree looking for a `potato_config.yaml` file.
Returns the project root path if found, NULL otherwise.

## Usage

``` r
find_potato_sack(path = NULL)
```

## Arguments

- path:

  Character. Starting path to search from. Defaults to current working
  directory.

## Value

Character path to the project root, or NULL if not inside a potato sack.
