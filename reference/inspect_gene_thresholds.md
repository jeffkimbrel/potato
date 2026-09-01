# Inspect gene-level annotation hits and thresholds

Shows all annotation hits for genes in a pathway, with scores vs
thresholds. Helps identify which genes are failing thresholds and by how
much.

## Usage

``` r
inspect_gene_thresholds(sack, potato_name = NULL, show_passing = FALSE)
```

## Arguments

- sack:

  PotatoSack object with annotation results

- potato_name:

  Optional: filter to specific pathway

- show_passing:

  Include hits that pass thresholds (default: FALSE)

## Value

Tibble with columns: potato, gene_id, gene_name, genome, tool, score,
threshold, margin, passed
