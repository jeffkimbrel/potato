#' Add message to PotatoSack
#'
#' Internal function to collect messages, warnings, and errors during workflow.
#' Messages are stored in the sack and can be retrieved with sack_messages().
#'
#' @param sack PotatoSack object
#' @param level Message level: "info", "warning", "error"
#' @param stage Workflow stage (e.g., "annotation", "scoring", "kofam")
#' @param message Message text
#' @param genome Optional genome name for per-genome messages
#' @param details Optional list of additional details
#'
#' @return Modified PotatoSack (invisibly)
#' @keywords internal
add_sack_message <- function(sack, level = c("info", "warning", "error"),
                             stage, message, genome = NULL, details = NULL) {
  level <- match.arg(level)

  # Trim whitespace from message
  message <- trimws(message)

  msg_record <- list(
    level = level,
    stage = stage,
    message = message,
    timestamp = Sys.time()
  )

  if (!is.null(genome)) {
    msg_record$genome <- genome
  }

  if (!is.null(details)) {
    msg_record$details <- details
  }

  # Append to messages list
  sack@messages <- c(sack@messages, list(msg_record))

  invisible(sack)
}


#' Print and collect message
#'
#' Helper function that both prints a message (if verbose) and collects it in the sack.
#'
#' @param sack PotatoSack object
#' @param message Message text
#' @param level Message level: "info", "warning", "error", "success", "header"
#' @param stage Workflow stage
#' @param genome Optional genome name
#' @param verbose Print to console (default TRUE)
#'
#' @return Modified PotatoSack (invisibly)
#' @keywords internal
sack_msg <- function(sack, message, level = c("info", "warning", "error", "success", "header"),
                     stage, genome = NULL, verbose = TRUE) {
  level <- match.arg(level)

  # Map display levels to storage levels
  storage_level <- switch(level,
    "header" = "info",
    "success" = "info",
    level
  )

  # Store message
  sack <- add_sack_message(sack, level = storage_level, stage = stage,
                          message = message, genome = genome)

  # Print if verbose
  if (verbose) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      # Plain text output
      if (level == "header") {
        cat("\n=== ", message, " ===\n", sep = "")
      } else if (level == "warning") {
        message("Warning: ", message)
      } else if (level == "error") {
        message("Error: ", message)
      } else {
        message(message)
      }
    } else {
      # CLI output
      if (level == "header") {
        cli::cli_h2(message)
      } else if (level == "warning") {
        cli::cli_alert_warning(message)
      } else if (level == "error") {
        cli::cli_alert_danger(message)
      } else if (level == "success") {
        cli::cli_alert_success(message)
      } else {
        cli::cli_alert_info(message)
      }
    }
  }

  invisible(sack)
}


#' Retrieve messages from PotatoSack
#'
#' Returns collected messages from the workflow. Useful for reviewing warnings
#' and errors that occurred during annotation or scoring.
#'
#' @param sack PotatoSack object
#' @param level Optional filter by message level: "info", "warning", "error"
#' @param stage Optional filter by workflow stage
#' @param genome Optional filter by genome name
#' @param as_dataframe Return as data frame instead of list (default TRUE)
#'
#' @return Data frame or list of messages
#' @export
#'
#' @examples
#' \dontrun{
#' # Get all messages
#' sack_messages(sack)
#'
#' # Get only warnings
#' sack_messages(sack, level = "warning")
#'
#' # Get messages from annotation stage
#' sack_messages(sack, stage = "annotation")
#'
#' # Get messages for specific genome
#' sack_messages(sack, genome = "genome001")
#' }
sack_messages <- function(sack, level = NULL, stage = NULL, genome = NULL,
                         as_dataframe = TRUE) {

  if (length(sack@messages) == 0) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("No messages collected")
    } else {
      cli::cli_alert_info("No messages collected")
    }
    return(if (as_dataframe) data.frame() else list())
  }

  msgs <- sack@messages

  # Filter by level
  if (!is.null(level)) {
    msgs <- Filter(function(m) m$level == level, msgs)
  }

  # Filter by stage
  if (!is.null(stage)) {
    msgs <- Filter(function(m) m$stage == stage, msgs)
  }

  # Filter by genome
  if (!is.null(genome)) {
    msgs <- Filter(function(m) !is.null(m$genome) && m$genome == genome, msgs)
  }

  if (length(msgs) == 0) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("No matching messages found")
    } else {
      cli::cli_alert_info("No matching messages")
    }
    return(if (as_dataframe) data.frame() else list())
  }

  if (!as_dataframe) {
    return(msgs)
  }

  # Convert to data frame
  df <- do.call(rbind, lapply(msgs, function(m) {
    data.frame(
      level = m$level,
      stage = m$stage,
      genome = if (!is.null(m$genome)) m$genome else NA_character_,
      message = m$message,
      timestamp = as.character(m$timestamp),
      stringsAsFactors = FALSE
    )
  }))

  if (!requireNamespace("tibble", quietly = TRUE)) {
    df
  } else {
    tibble::as_tibble(df)
  }
}


#' Print summary of sack messages
#'
#' Displays a formatted summary of messages collected during workflow.
#'
#' @param sack PotatoSack object
#' @param max_show Maximum number of messages to show per category (default 10)
#'
#' @return Invisibly returns the sack
#' @export
#'
#' @examples
#' \dontrun{
#' sack_message_summary(sack)
#' }
sack_message_summary <- function(sack, max_show = 10) {

  if (length(sack@messages) == 0) {
    if (!requireNamespace("cli", quietly = TRUE)) {
      message("No messages collected")
    } else {
      cli::cli_alert_info("No messages collected")
    }
    return(invisible(sack))
  }

  # Count by level
  levels <- sapply(sack@messages, function(m) m$level)
  level_counts <- table(levels)

  if (!requireNamespace("cli", quietly = TRUE)) {
    cat("\n=== Message Summary ===\n")
    cat("Total messages:", length(sack@messages), "\n")
    for (lvl in names(level_counts)) {
      cat("  ", lvl, ": ", level_counts[lvl], "\n", sep = "")
    }
  } else {
    cli::cli_h2("Message Summary")
    cli::cli_text("Total: {length(sack@messages)}")
    cli::cli_ul()
    for (lvl in names(level_counts)) {
      symbol <- switch(lvl,
        "error" = cli::symbol$cross,
        "warning" = cli::symbol$warning,
        "info" = cli::symbol$info,
        cli::symbol$bullet
      )
      cli::cli_li("{symbol} {lvl}: {level_counts[lvl]}")
    }
    cli::cli_end()
  }

  # Show recent messages by level
  for (lvl in c("error", "warning", "info")) {
    msgs <- Filter(function(m) m$level == lvl, sack@messages)
    if (length(msgs) == 0) next

    n_show <- min(length(msgs), max_show)
    msgs_show <- utils::tail(msgs, n_show)

    if (!requireNamespace("cli", quietly = TRUE)) {
      cat("\n", toupper(lvl), "S:\n", sep = "")
      for (msg in msgs_show) {
        cat("  [", msg$stage, "] ", msg$message, "\n", sep = "")
        if (!is.null(msg$genome)) {
          cat("    Genome: ", msg$genome, "\n", sep = "")
        }
      }
      if (length(msgs) > max_show) {
        cat("  ... and ", length(msgs) - max_show, " more\n", sep = "")
      }
    } else {
      cli::cli_h3(paste0(toupper(lvl), "S"))
      for (msg in msgs_show) {
        genome_txt <- if (!is.null(msg$genome)) paste0(" [", msg$genome, "]") else ""
        cli::cli_text("{.field {msg$stage}}{genome_txt}: {msg$message}")
      }
      if (length(msgs) > max_show) {
        cli::cli_text("{cli::symbol$ellipsis} and {length(msgs) - max_show} more")
      }
    }
  }

  cat("\n")
  if (!requireNamespace("cli", quietly = TRUE)) {
    cat("Use sack_messages(sack) to view all messages as a data frame\n")
  } else {
    cli::cli_alert_info("Use {.code sack_messages(sack)} to view all messages")
  }

  invisible(sack)
}
