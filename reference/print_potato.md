# Print potato structure as text

Shows pathway flow with genes, compounds, and structure

## Usage

``` r
print_potato(
  potato,
  compact = TRUE,
  show_compounds = TRUE,
  show_databases = FALSE,
  show_ec = FALSE
)
```

## Arguments

- potato:

  Potato S7 object, v2 potato object, or path to JSON

- compact:

  Show compact one-line view (default: TRUE)

- show_compounds:

  Include compound names in flow (default: TRUE, ignored if compact)

- show_databases:

  Show database annotations (kofam, blast, hmm, etc.) (default: FALSE)

- show_ec:

  Show EC numbers (default: FALSE)
