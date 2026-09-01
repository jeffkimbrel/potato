# Align node coordinates to grid or by axis

Snap coordinates to a grid or align nodes that are nearly
horizontal/vertical. Useful for cleaning up manually-arranged layouts.

## Usage

``` r
align_coordinates(
  coords_path,
  output_path = NULL,
  grid_size = NULL,
  align_threshold = 2.5,
  axis = "both"
)
```

## Arguments

- coords_path:

  Path to coordinates JSON file

- output_path:

  Path to save aligned coordinates (default: overwrites input)

- grid_size:

  Snap to grid of this size (NULL = no grid snapping)

- align_threshold:

  Percentage of coordinate range for alignment (NULL = no alignment).
  Nodes within this percentage distance are aligned to same axis.
  Examples: 2.5 = gentle (default), 5.0 = moderate, 10.0 = aggressive.
  Scale-independent - works regardless of coordinate system.

- axis:

  Axis to align: "y" (horizontal lines), "x" (vertical lines), or "both"
