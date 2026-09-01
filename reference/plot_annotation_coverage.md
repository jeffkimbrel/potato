# Plot annotation coverage heatmap

Shows which pathways have been annotated with which databases (kofam,
blast, hmm). Helps identify gaps in annotation coverage.

## Usage

``` r
plot_annotation_coverage(sack)
```

## Arguments

- sack:

  PotatoSack object with provenance data

## Value

ggplot2 heatmap

## Details

Cell colors:

- Green: Pathway was checked with this database (pathway has genes AND
  was in potatoes_requested)

- Gray: Pathway has no genes for this database (not applicable)

- Red: Pathway HAS genes for this database but was NOT checked
  (annotation gap!)
