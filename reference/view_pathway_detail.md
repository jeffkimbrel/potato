# View pathway details as HTML table in browser

Opens a simple HTML table in the system browser showing pathway gene
details. Useful when terminal output is too wide to read comfortably.

## Usage

``` r
view_pathway_detail(potato, pathway = NULL, layout = "fr")
```

## Arguments

- potato:

  Potato S7 object or path to JSON

- pathway:

  For multi-pathway networks, which pathway to show (single pathway ID
  only)

- layout:

  Layout algorithm for visualization: "xy" (curated coords), "fr", "kk",
  "sugiyama", "tree", "circle", "grid" (default: "fr")
