# Get simple pathway verification status table

Returns a clean tibble with one row per pathway showing verification
status.

## Usage

``` r
get_pathway_verification(input)
```

## Arguments

- input:

  Either a PotatoSack object or a directory path containing potato JSON
  files

## Value

Tibble with columns: potato, pathway, verified
