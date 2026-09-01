# Run kofam annotation on all genomes

Run kofam annotation on all genomes

## Usage

``` r
run_kofam(
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

  If FALSE (default), error if kofam results already exist

## Value

Modified PotatoSack with kofam results in @results\$kofam
