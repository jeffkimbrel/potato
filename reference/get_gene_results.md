# Get gene-level annotation results

Returns a tibble of all gene-level annotation hits across all tools,
formatted as a long table compatible with GATOR output format. Includes
a "passed" column indicating whether each hit passes quality thresholds.

## Usage

``` r
get_gene_results(sack)
```

## Arguments

- sack:

  PotatoSack object with annotation results

## Value

Tibble with columns: genome, potato, node_id, gene_id, tool, passed, and
tool-specific columns
