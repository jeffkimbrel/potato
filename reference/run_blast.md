# Run BLAST annotation on all genomes

Run BLAST annotation on all genomes

## Usage

``` r
run_blast(
  sack,
  potato_names = NULL,
  conda_env = NULL,
  workers = NULL,
  overwrite = FALSE
)
```

## Arguments

- sack:

  PotatoSack object

- potato_names:

  Vector of potato names (NULL = all)

- conda_env:

  Optional conda environment name (defaults to config setting)

- workers:

  Number of parallel workers (defaults to config setting, 1 =
  sequential)

- overwrite:

  If FALSE (default), error if blast results already exist

## Value

Modified PotatoSack with blast results in @results\$blast
