# Get detailed blast annotation results

Returns detailed information about blast annotation hits, including
per-genome and per-potato summaries and plots.

## Usage

``` r
get_blast_details(sack)
```

## Arguments

- sack:

  PotatoSack object

## Value

A list containing:

- summary:

  Tibble with per-genome hit counts

- per_potato:

  Tibble with per-potato hit statistics

- plot:

  ggplot object showing hits per genome

## Examples

``` r
if (FALSE) { # \dontrun{
details <- get_blast_details(sack)
details$summary
details$per_potato
details$plot
} # }
```
