# Get pathway-level scoring results

Returns a tibble of pathway completion scores across all genomes,
formatted for analysis and export. Includes potato_hash for version
tracking.

## Usage

``` r
get_pathway_scores(sack)
```

## Arguments

- sack:

  PotatoSack object with scoring results

## Value

Tibble with columns: genome, potato, potato_hash, pathway_name,
total_steps_detected, total_steps, fraction, min_fraction, present,
essential_total_steps_detected, essential_steps, essential_fraction,
essential_pathway_present
