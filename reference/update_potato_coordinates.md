# Update potato JSON with node coordinates

After arranging nodes in visNetwork and exporting coordinates, use this
function to add x,y fields to your potato JSON.

## Usage

``` r
update_potato_coordinates(
  potato_path,
  coords_path,
  output_path = NULL,
  with_compounds = NULL
)
```

## Arguments

- potato_path:

  Path to potato JSON file

- coords_path:

  Path to coordinates JSON file (from visNetwork export)

- output_path:

  Path to save updated potato JSON (default: overwrites original)

- with_compounds:

  Logical - if NULL (default), auto-detects from coordinate file; if
  TRUE/FALSE, forces that mode
