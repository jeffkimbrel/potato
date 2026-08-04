#' Align node coordinates to grid or by axis
#'
#' Snap coordinates to a grid or align nodes that are nearly horizontal/vertical.
#' Useful for cleaning up manually-arranged layouts.
#'
#' @param coords_path Path to coordinates JSON file
#' @param output_path Path to save aligned coordinates (default: overwrites input)
#' @param grid_size Snap to grid of this size (NULL = no grid snapping)
#' @param align_threshold Nodes within this distance are aligned to same axis (NULL = no alignment)
#' @param axis Axis to align: "y" (horizontal lines), "x" (vertical lines), or "both"
#'
#' @export
align_coordinates <- function(coords_path, output_path = NULL,
                              grid_size = NULL, align_threshold = 10, axis = "both") {

  if (is.null(output_path)) {
    output_path <- coords_path
  }

  # Read coordinates
  coords <- jsonlite::read_json(coords_path, simplifyVector = TRUE)

  original_coords <- coords

  # Grid snapping
  if (!is.null(grid_size)) {
    coords$x <- round(coords$x / grid_size) * grid_size
    coords$y <- round(coords$y / grid_size) * grid_size
    cli::cli_alert_info("Snapped to {grid_size}-unit grid")
  }

  # Axis alignment
  if (!is.null(align_threshold)) {

    # Align Y (horizontal lines)
    if (axis %in% c("y", "both")) {
      coords_sorted <- coords[order(coords$y), ]
      groups <- list()
      current_group <- 1
      groups[[current_group]] <- 1

      for (i in 2:nrow(coords_sorted)) {
        if (abs(coords_sorted$y[i] - coords_sorted$y[i-1]) <= align_threshold) {
          groups[[current_group]] <- c(groups[[current_group]], i)
        } else {
          current_group <- current_group + 1
          groups[[current_group]] <- i
        }
      }

      # Align each group to median Y
      y_changes <- 0
      for (group in groups) {
        if (length(group) > 1) {
          group_y <- coords_sorted$y[group]
          median_y <- median(group_y)
          coords_sorted$y[group] <- median_y
          y_changes <- y_changes + 1

          node_ids <- coords_sorted$id[group]
          cli::cli_alert_success("Aligned {length(group)} nodes horizontally at Y={median_y}: {paste(node_ids, collapse=', ')}")
        }
      }

      # Restore original order
      coords <- coords_sorted[match(coords$id, coords_sorted$id), ]
    }

    # Align X (vertical lines)
    if (axis %in% c("x", "both")) {
      coords_sorted <- coords[order(coords$x), ]
      groups <- list()
      current_group <- 1
      groups[[current_group]] <- 1

      for (i in 2:nrow(coords_sorted)) {
        if (abs(coords_sorted$x[i] - coords_sorted$x[i-1]) <= align_threshold) {
          groups[[current_group]] <- c(groups[[current_group]], i)
        } else {
          current_group <- current_group + 1
          groups[[current_group]] <- i
        }
      }

      # Align each group to median X
      x_changes <- 0
      for (group in groups) {
        if (length(group) > 1) {
          group_x <- coords_sorted$x[group]
          median_x <- median(group_x)
          coords_sorted$x[group] <- median_x
          x_changes <- x_changes + 1

          node_ids <- coords_sorted$id[group]
          cli::cli_alert_success("Aligned {length(group)} nodes vertically at X={median_x}: {paste(node_ids, collapse=', ')}")
        }
      }

      # Restore original order
      coords <- coords_sorted[match(coords$id, coords_sorted$id), ]
    }
  }

  # Write aligned coordinates
  jsonlite::write_json(coords, output_path, pretty = TRUE, auto_unbox = TRUE)

  # Report changes
  x_diff <- sum(abs(coords$x - original_coords$x) > 0.01)
  y_diff <- sum(abs(coords$y - original_coords$y) > 0.01)

  cli::cli_alert_success("Aligned coordinates saved to: {.file {output_path}}")
  cli::cli_alert_info("Changed: {x_diff} X coordinates, {y_diff} Y coordinates")

  invisible(output_path)
}
