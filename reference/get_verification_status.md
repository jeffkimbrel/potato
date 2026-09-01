# Get verification status of pathways

Shows which pathways in potatoes are verified or need verification. Can
check a single potato or all potatoes in a sack. Optionally blocks
editing of verified pathways (guard mode).

## Usage

``` r
get_verification_status(
  potato_or_sack,
  pathway_id = NULL,
  abort_if_verified = FALSE,
  force = FALSE
)
```

## Arguments

- potato_or_sack:

  Potato object, potato path, PotatoSack object, or list of potatoes

- pathway_id:

  Optional pathway ID to check within a potato

- abort_if_verified:

  Logical. If TRUE, aborts with error if pathway is verified (guard
  mode). Default FALSE.

- force:

  Logical. If TRUE, bypasses verification guard (only used when
  abort_if_verified=TRUE). Default FALSE.

## Value

Prints verification summary and invisibly returns a tibble (or TRUE in
guard mode)
