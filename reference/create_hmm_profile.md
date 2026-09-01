# Create filtered HMM profile file from potato detection terms

Extracts only the HMM profiles needed by potatoes from the configured
HMM files and creates a filtered HMM profile file (may be concatenated).
Also extracts trusted cutoff (TC) values from profiles for use in
scoring.

## Usage

``` r
create_hmm_profile(sack, potato_names = NULL)
```

## Arguments

- sack:

  PotatoSack object

- potato_names:

  Vector of potato names (NULL = all)

## Value

List with hmm_profile path, tc_values, and modified sack
