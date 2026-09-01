# Summarize hmm annotation

Displays messages from run_hmm(). Messages are retrieved from provenance
data stored in the sack, so this works even if the original run_hmm()
call was in an eval=FALSE chunk.

## Usage

``` r
summarize_hmm(sack, verbose = FALSE)
```

## Arguments

- sack:

  PotatoSack object

- verbose:

  Logical. If TRUE, prints detailed messages (default: TRUE)

## Value

A list containing:

- summary:

  Tibble with annotation summary (n_genomes, n_potatoes, n_profiles,
  tool_version)

- messages:

  Tibble with columns type, message

- status:

  List with ok (logical)

## Examples

``` r
if (FALSE) { # \dontrun{
sack <- run_hmm(sack)
result <- summarize_hmm(sack)
result$summary
result$messages  # Render this tibble in qmd
} # }
```
