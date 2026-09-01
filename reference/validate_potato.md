# Validate potato structure

Validates a potato JSON structure for common errors and required fields.
Handles both multi-pathway networks (new schema) and single-pathway
potatoes (legacy schema). Users should run this on custom potatoes
before using them for annotation.

## Usage

``` r
validate_potato(potato, strict = FALSE)
```

## Arguments

- potato:

  Potato object or list (raw JSON data)

- strict:

  Logical. If TRUE, performs additional strict checks (default FALSE)

## Value

List with 'valid' (logical), 'errors' (character vector), and 'warnings'
(character vector)
