# Load potato configuration from YAML

Reads and validates a potato_config.yaml file. Returns a structured list
with tool configurations, database paths, and thresholds.

## Usage

``` r
load_potato_config(config_path = NULL)
```

## Arguments

- config_path:

  Path to potato_config.yaml file. If NULL, searches for it using
  find_potato_sack() from current directory.

## Value

List with parsed configuration

## Details

The configuration file defines:

- Tool executables and database paths (kofam, pfam, hmmer, blast)

- Default thresholds for each tool

- Project paths and settings

## Examples

``` r
if (FALSE) { # \dontrun{
config <- load_potato_config("potato_config.yaml")
config <- load_potato_config()  # Auto-finds in current sack
} # }
```
