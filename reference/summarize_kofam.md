# Summarize kofam annotation

Displays messages from run_kofam(). Messages are retrieved from
provenance data stored in the sack, so this works even if the original
run_kofam() call was in an eval=FALSE chunk.

## Usage

``` r
summarize_kofam(sack, verbose = FALSE)
```

## Arguments

- sack:

  PotatoSack object

- verbose:

  Logical. If TRUE, prints detailed messages (default: TRUE)

## Value

A list containing:

- summary:

  Tibble with annotation summary (n_genomes, n_potatoes, n_kos,
  tool_version)

- messages:

  Tibble with columns type, message

- status:

  List with ok (logical)

## Examples

``` r
if (FALSE) { # \dontrun{
sack <- run_kofam(sack)
result <- summarize_kofam(sack)
result$summary
result$messages  # Render this tibble in qmd
} # }
```
