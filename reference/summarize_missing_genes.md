# Identify genes missing across genomes

Analyzes which genes are systematically missing, with distinction
between genes with no hits, hits below threshold, and hits passing
threshold.

## Usage

``` r
summarize_missing_genes(sack, potato_name = NULL, min_genomes = 3)
```

## Arguments

- sack:

  PotatoSack object with annotation results and scores

- potato_name:

  Optional: filter to specific pathway

- min_genomes:

  Minimum genomes to consider (default: 3)

## Value

Tibble showing genes and detection status across thresholds
