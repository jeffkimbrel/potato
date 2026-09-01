# Find pathways that are close to being called present

Identifies "near miss" pathways that are just below the threshold, which
may indicate threshold tuning issues.

## Usage

``` r
find_near_miss_pathways(sack, buffer = 0.1)
```

## Arguments

- sack:

  PotatoSack object with scores

- buffer:

  Distance from threshold to consider "near" (default: 0.1)

## Value

Tibble of near-miss pathways
