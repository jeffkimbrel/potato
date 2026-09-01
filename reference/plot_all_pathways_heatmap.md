# Plot comprehensive pathway completion heatmap across all genomes

Creates a heatmap showing pathway completion for all potatoes and
genomes. Tiles show fraction of pathway detected, optionally normalized
by threshold.

## Usage

``` r
plot_all_pathways_heatmap(
  sack,
  normalize_by_threshold = TRUE,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_labels = TRUE,
  clustering_method = "complete"
)
```

## Arguments

- sack:

  PotatoSack object with scores

- normalize_by_threshold:

  Logical. If TRUE, shows fraction/min_fraction ratio with color
  breakpoint at 1.0 (detection threshold). If FALSE, shows raw fraction.

- cluster_rows:

  Logical. Cluster pathways by similarity (default: FALSE)

- cluster_cols:

  Logical. Cluster genomes by similarity (default: TRUE)

- show_labels:

  Logical. Show pathway and genome labels (default: TRUE)

- clustering_method:

  Clustering method for hclust: "complete", "average", "single",
  "ward.D2" (default: "complete")

## Value

A ggplot2 object
