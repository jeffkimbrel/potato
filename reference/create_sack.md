# Create a PotatoSack from project directory

Constructs a PotatoSack S7 object from a project directory containing
potato_config.yaml and potatoes. To save/load sacks, use standard R
functions: saveRDS(sack, "sack.rds") and readRDS("sack.rds").

## Usage

``` r
create_sack(path = NULL)
```

## Arguments

- path:

  Character. Path to sack directory. If NULL, searches upward from
  current directory.

## Value

PotatoSack S7 object
