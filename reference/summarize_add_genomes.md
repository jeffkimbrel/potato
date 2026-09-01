# Summarize genome addition history

Displays messages from add_genomes() calls. Messages are reconstructed
from provenance data stored in the sack, so this works even if the
original add_genomes() call was in an eval=FALSE chunk.

## Usage

``` r
summarize_add_genomes(sack, verbose = FALSE)
```

## Arguments

- sack:

  PotatoSack object

- verbose:

  Logical. If TRUE, prints detailed messages (default: TRUE)

## Value

A list containing:

- summary:

  Tibble with per-call statistics (call_number, timestamp, n_added,
  total_after)

- messages:

  Tibble with columns type, message

- status:

  List with ok (logical), total_genomes (integer)

## Examples

``` r
if (FALSE) { # \dontrun{
sack <- add_genomes(sack, "~/genomes/*.faa")
result <- summarize_add_genomes(sack)
result$summary
result$messages  # Render this tibble in qmd
} # }
```
