# Compute MD5 hash of potato R object

Computes a hash of the potato's functional content (nodes, edges,
scoring, etc.) to track which version was used for annotation. Excludes
metadata fields like json_path and graph cache that don't affect
annotation behavior.

## Usage

``` r
compute_potato_hash(potato)
```

## Arguments

- potato:

  Potato S7 object

## Value

Character string with MD5 hash

## Details

Specifically excludes cosmetic fields that don't affect
detection/scoring:

- x, y, x_compounds, y_compounds (visualization coordinates)

- notes (documentation only)
