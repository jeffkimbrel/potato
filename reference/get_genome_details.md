# Get detailed genome information

Returns a list with detailed information about all genomes in the sack,
including file sizes, protein counts, and QC plots.

## Usage

``` r
get_genome_details(sack)
```

## Arguments

- sack:

  PotatoSack object

## Value

A list containing:

- summary:

  Tibble with genome information (genome, file_path, file_type, md5,
  file_size_mb, n_proteins, added_in_call)

- plot:

  ggplot object with protein count bar plot

## Examples

``` r
if (FALSE) { # \dontrun{
result <- get_genome_details(sack)
result$summary
result$plot
} # }
```
