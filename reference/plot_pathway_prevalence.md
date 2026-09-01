# Plot pathway prevalence across genomes

Creates a horizontal bar chart showing how many genomes have each
pathway. Useful for identifying widespread vs rare metabolic
capabilities.

## Usage

``` r
plot_pathway_prevalence(
  sack,
  min_genomes = 0,
  pathway_type = NULL,
  sort_by = "count"
)
```

## Arguments

- sack:

  PotatoSack object with scoring results

- min_genomes:

  Minimum number of genomes to include pathway (default: 0)

- pathway_type:

  For multi-pathway networks, filter to "variant" or "independent"
  (default: NULL = all)

- sort_by:

  Sort bars by "count" (default), "name", or "fraction"

## Value

A ggplot2 object
