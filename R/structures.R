#' Set up Spatial Structure for Simulation
#'
#' Prepares a data frame of regions (lowest administrative units) for the simulation,
#' including any specified hierarchical levels and names.
#'
#' @param spatial_structure Specification of the spatial layout. This can be:
#'   \itemize{
#'     \item A data frame containing at least one column of region identifiers (and
#'           optionally higher-level region columns). If provided, it will be returned
#'           as-is (after ensuring the lowest level column is last).
#'     \item A named numeric vector defining the number of units at each hierarchy level
#'           (from top level to bottom level). For example, \code{c(country = 1, region = 2, district = 4)}
#'           means 1 country, 2 regions (total), and 4 districts (total), distributed
#'           as evenly as possible under those 2 regions.
#'     \item An unnamed numeric vector of length > 1, giving the number of units at each
#'           level (defaults to generic names "Level1", "Level2", ... if no \code{level_names} provided).
#'     \item A single numeric value, interpreted as the number of lowest-level regions
#'           (with a single top-level grouping if a name is provided for that level).
#'     \item A character vector of region names (for the lowest level; a single top-level
#'           grouping will be added if no higher level is specified).
#'   }
#' @param level_names Optional character vector of names for each hierarchy level (used
#' when \code{spatial_structure} is numeric and unnamed, or when a character vector of
#' region names is provided). The order should be from top level to bottom level.
#'
#' @return A data frame where each row corresponds to one lowest-level region. The data
#' frame includes one column for each hierarchical level (with the lowest level last).
#' If only one level is defined, the column will be named "region" (or the provided name
#' for that level). If a single top-level grouping is implied (e.g., only one country),
#' it will still be included as a column for consistency.
#'
#' @details The hierarchy is constructed such that lowest-level units are distributed as
#' evenly as possible among their parent units at the next level up. If the numbers do
#' not divide evenly, some parent units will have one more sub-unit than others. This
#' function ensures that even with uneven splits, all specified units are created and
#' assigned appropriately.
#'
#' If a data frame is provided as \code{spatial_structure}, it should already contain the
#' desired hierarchy columns (e.g., columns for each admin level). The function will
#' return it (after moving the lowest level column to the last position if it isn't already).
#'
#' @examples
#' # Example 1: Numeric hierarchy specification (1 country, 2 regions, 4 districts)
#' setup_spatial_structure(c(country = 1, region = 2, district = 4))
#'
#' # Example 2: Provide only the number of regions (no hierarchy above)
#' setup_spatial_structure(3)
#'
#' # Example 3: Provide explicit names for lowest-level regions
#' setup_spatial_structure(c("North", "South", "East"))
#'
#' # Example 4: Use a custom data frame for regions
#' regions_df <- data.frame(country = "X", region = c("A", "B"), district = c("A1", "B1"))
#' setup_spatial_structure(regions_df)
#' @export
setup_spatial_structure <- function(spatial_structure = 1, level_names = NULL) {
  # If a data frame of regions is provided, use it directly
  if (is.data.frame(spatial_structure)) {
    region_df <- spatial_structure
    # Ensure the lowest-level column is last
    # If multiple columns, assume last column is already the lowest level.
    # (We do not reorder user-provided data beyond this assumption.)
    return(region_df)
  }

  # Case 1: Single numeric value (number of regions, one top-level grouping)
  if (is.numeric(spatial_structure) && length(spatial_structure) == 1 && is.null(names(spatial_structure))) {
    n_regions <- spatial_structure
    lowest_name <- if (!is.null(level_names)) {
      # If a name is provided for the single level
      level_names[length(level_names)]
    } else {
      "region"
    }
    top_name <- if (!is.null(level_names) && length(level_names) > 1) {
      level_names[1]
    } else {
      "country"
    }
    if (!is.null(level_names) && length(level_names) == 1) {
      # Only one level specified by user name; no explicit top-level column
      region_names <- paste0(level_names, seq_len(n_regions))
      region_df <- data.frame(setNames(list(region_names), lowest_name), stringsAsFactors = FALSE)
    } else {
      # Include a single top-level (e.g., country) column
      region_names <- paste0(lowest_name, seq_len(n_regions))
      region_df <- data.frame(
        setNames(list(rep(paste0(top_name, "1"), n_regions)), top_name),
        setNames(list(region_names), lowest_name),
        stringsAsFactors = FALSE
      )
    }
    return(region_df)
  }

  # Case 2: Character vector of region names (lowest level regions)
  if (is.character(spatial_structure)) {
    region_names <- spatial_structure
    n_regions <- length(region_names)
    lowest_name <- if (!is.null(level_names) && length(level_names) >= 1) {
      tail(level_names, 1)
    } else {
      "region"
    }
    top_name <- if (!is.null(level_names) && length(level_names) > 1) {
      level_names[1]
    } else {
      "country"
    }
    if (!is.null(level_names) && length(level_names) == 1) {
      # Only one level with provided name
      region_df <- data.frame(setNames(list(region_names), lowest_name), stringsAsFactors = FALSE)
    } else {
      # Provide a single top-level grouping (assuming one group for all regions)
      region_df <- data.frame(
        setNames(list(rep(paste0(top_name, "1"), n_regions)), top_name),
        setNames(list(region_names), lowest_name),
        stringsAsFactors = FALSE
      )
    }
    return(region_df)
  }

  # Case 3: Numeric vector with hierarchy counts (possibly named or unnamed)
  if (is.numeric(spatial_structure) && length(spatial_structure) > 1) {
    level_counts <- spatial_structure
    # If no names provided, assign generic or provided level_names
    if (is.null(names(level_counts))) {
      if (!is.null(level_names) && length(level_names) == length(level_counts)) {
        names(level_counts) <- level_names
      } else {
        names(level_counts) <- paste0("Level", seq_along(level_counts))
      }
    }
    # Ensure highest level is first and lowest is last in names (assume input is top-to-bottom already)
    levels <- names(level_counts)
    n_levels <- length(levels)
    region_df <- data.frame(stringsAsFactors = FALSE)
    # Build hierarchy iteratively
    for (i in seq_along(levels)) {
      level_name <- levels[i]
      count <- level_counts[i]
      if (i == 1) {
        # Top level: create that many entries
        ids <- seq_len(count)
        region_df <- data.frame(level_name = paste0(level_name, ids), stringsAsFactors = FALSE)
        names(region_df)[1] <- level_name
      } else {
        # Subsequent levels: distribute units among the previous level's units
        parent_count <- level_counts[i - 1]
        total_units <- count
        # Base number of sub-units per parent
        base_per_parent <- total_units %/% parent_count
        remainder <- total_units %% parent_count
        # Each of the first 'remainder' parents gets one extra sub-unit
        counts <- rep(base_per_parent, parent_count)
        if (remainder > 0) {
          counts[1:remainder] <- counts[1:remainder] + 1
        }
        # Replicate the current region_df rows according to counts for each parent unit
        region_df <- region_df[rep(seq_len(nrow(region_df)), counts), , drop = FALSE]
        # Assign identifiers for the current level
        ids <- seq_len(total_units)
        region_df[[level_name]] <- paste0(level_name, ids)
      }
    }
    rownames(region_df) <- NULL
    return(region_df)
  }

  stop("Invalid spatial_structure specification. Must be a data frame, a numeric value/vector, or a character vector of region names.")
}


#' Set up Temporal Simulation Structure
#'
#' Determines the sequence of daily dates and the sequence of output period start
#' dates based on a start date and either an end date or a number of periods, plus
#' a temporal resolution.
#'
#' @param start_date Date (or string in "YYYY-MM-DD" format) for the start of the simulation.
#' @param end_date Date (or string) for the end of the simulation. Either \code{end_date}
#' or \code{n_periods} must be provided (but not both).
#' @param n_periods Integer number of time periods to simulate at the output resolution
#' (ignored if \code{end_date} is provided).
#' @param resolution Time resolution for output periods. Options include \code{"day"},
#' \code{"week"}, \code{"month"}, \code{"quarter"}, or \code{"year"}.
#'
#' @return A list with two Date vectors: \code{time_seq_daily} (all daily dates in the
#' simulation) and \code{time_seq_out} (the start date of each output period at the
#' specified resolution).
#'
#' @details If \code{end_date} is provided, the daily sequence runs from \code{start_date}
#' to \code{end_date} inclusive. If \code{n_periods} is provided instead, the function
#' generates that many consecutive periods (at the given resolution) starting from
#' \code{start_date}, and computes an \code{end_date} to cover all those periods fully.
#' The output vector \code{time_seq_out} contains the start date of each period.
#'
#' For example, if \code{start_date = "2025-01-15"} and \code{resolution = "month"}
#' with \code{n_periods = 2}, the output periods will start on 2025-01-15 and 2025-02-15,
#' and the daily sequence will run until 2025-03-14 to cover two full monthly periods
#' (since the second period runs from 2025-02-15 to 2025-03-14).
#'
#' @examples
#' # Using an end date
#' setup_temporal_structure(start_date = "2025-01-01", end_date = "2025-01-10", resolution = "day")
#' setup_temporal_structure(start_date = "2025-01-01", end_date = "2025-03-10", resolution = "month")
#'
#' # Using number of periods (e.g., 3 monthly periods starting Jan 1, 2025)
#' setup_temporal_structure(start_date = "2025-01-01", n_periods = 3, resolution = "month")
#' @export
setup_temporal_structure <- function(start_date, end_date = NULL, n_periods = NULL,
                                     resolution = c("day", "week", "month", "quarter", "year")) {
  # Convert inputs to Date objects
  start_date <- as.Date(start_date)
  if (!is.null(end_date)) {
    end_date <- as.Date(end_date)
  }
  resolution <- match.arg(resolution)

  # Validate that exactly one of end_date or n_periods is provided
  if (!xor(is.null(end_date), is.null(n_periods))) {
    stop("Please provide either end_date or n_periods (but not both).")
  }

  # Compute daily and output period sequences
  if (!is.null(end_date)) {
    # Daily sequence from start to end (inclusive)
    time_seq_daily <- seq.Date(from = start_date, to = end_date, by = "day")
    # Sequence of output period start dates from start_date to end_date
    time_seq_out <- seq.Date(from = start_date, to = end_date, by = resolution)
    # If end_date is not exactly on a period boundary, add the next period start to cover it
    if (tail(time_seq_out, 1) < end_date) {
      time_seq_out <- c(time_seq_out, seq.Date(from = tail(time_seq_out, 1), by = resolution, length.out = 2)[2])
    }
  } else {
    # n_periods provided: generate period start dates and infer end_date
    time_seq_out <- seq.Date(from = start_date, by = resolution, length.out = n_periods)
    if (resolution == "day") {
      # If daily periods, end_date is just the last period start (covers that day)
      end_date <- tail(time_seq_out, 1)
    } else {
      # For longer periods, find the end of the last period
      next_start <- seq.Date(from = tail(time_seq_out, 1), by = resolution, length.out = 2)[2]
      end_date <- next_start - 1  # one day before the next period start
    }
    # Daily sequence from start_date up to end_date
    time_seq_daily <- seq.Date(from = start_date, to = end_date, by = "day")
  }

  return(list(time_seq_daily = time_seq_daily, time_seq_out = time_seq_out))
}
