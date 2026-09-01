# Generate KEGG URL for a specific pathway

Creates a KEGG database URL with all unique KO IDs and compound IDs from
a specific pathway. Useful for viewing pathway entries together.

## Usage

``` r
get_kegg_url(potato, pathway_id = NULL, genes = TRUE, compounds = TRUE)
```

## Arguments

- potato:

  Potato object or path to potato JSON

- pathway_id:

  Pathway ID within the potato (required for multi-pathway networks)

- genes:

  Logical. Include KO IDs from genes (default: TRUE)

- compounds:

  Logical. Include compound IDs (default: TRUE)

## Value

Character string with KEGG URL, or NULL if no IDs found
