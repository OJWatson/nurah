#' Introduce missing data (gaps) into a simulated dataset
#' @description
#' Apply specified missingness patterns to a dataset to simulate data gaps. Supports three types of missingness: MCAR (Missing Completely At Random), MAR (Missing At Random, via systematic gaps in certain regions or periods), and MNAR (Missing Not At Random, dependent on the data values). This function allows defining different missingness patterns for different variables.
#'
#' @param data A data frame containing the complete simulated data (e.g., output from `simulate_crisis_data`), with columns for time, spatial identifiers, and variables.
#' @param patterns A list of missingness pattern specifications. Each element of the list should be a list (or named list) with components describing one pattern. Required fields:
#'   \itemize{
#'     \item \code{variable}: Name of the variable (column in \code{data}) to introduce missingness in.
#'     \item \code{type}: Type of missingness, one of \code{"MCAR"}, \code{"MAR"}, or \code{"MNAR"}.
#'   }
#'   Depending on \code{type}, additional fields should be provided:
#'   \itemize{
#'     \item \strong{MCAR:} Provide \code{prop} (proportion of values to set missing at random). Optionally, can include \code{region} or \code{regions} and/or \code{start}, \code{end} to restrict random missingness to a specific spatial or temporal subset (otherwise missingness is applied randomly over all observations of the variable).
#'     \item \strong{MAR:} Provide criteria for systematic missingness. Typically \code{region} (a region name or vector of regions) and/or a time window via \code{start} and \code{end} (dates or strings). All values of the specified variable in the given region(s) and time range will be set to NA. If only \code{region} is given, the variable is entirely missing for those region(s). If only a time window is given, the variable is missing for that period in all regions.
#'     \item \strong{MNAR:} Provide \code{threshold} (numeric) and \code{direction} (">" or "<") to indicate that values above or below a threshold should be made missing. For example, \code{threshold=5, direction=">"} means any value greater than 5 in the specified variable will be replaced with NA (simulating not-at-random missingness dependent on the value itself). Optionally, \code{region} and/or time \code{start}/\code{end} can also be specified to apply this rule only in certain subsets.
#'   }
#'
#' @details
#' The function iterates over each pattern in \code{patterns} and introduces missing values accordingly. Multiple patterns can affect the same variable sequentially (e.g., a variable could have some baseline random missingness and an additional block of MAR missingness). Time values in \code{start} and \code{end} can be provided as Date objects or strings in "YYYY-MM-DD" format (they will be converted to Date). The \code{data} is not modified in place; a new data frame with the specified gaps is returned.
#'
#' @return A data frame identical to \code{data} but with specified entries replaced by \code{NA} to represent missing values.
#'
#' @examples
#' # Introduce 10% random missingness to malnutrition_rate and a systematic gap for mortality_rate
#' patterns <- list(
#'   list(variable="malnutrition_rate", type="MCAR", prop=0.1),
#'   list(variable="mortality_rate", type="MAR", region="RegionB", start="2020-06-01", end="2020-08-31")
#' )
#' sim_with_gaps <- introduce_data_gaps(sim_data, patterns)
#' summary(sim_with_gaps$malnutrition_rate)  # should show some NAs
#' # All mortality_rate values in RegionB (admin1) during Jun-Aug 2020 are NA:
#' all(is.na(subset(sim_with_gaps, admin1=="RegionB" & time >= as.Date("2020-06-01") & time <= as.Date("2020-08-31"))$mortality_rate))
#'
#' @export
introduce_data_gaps <- function(data, patterns) {
  data_out <- data
  for (pat in patterns) {
    if (is.null(pat$variable) || is.null(pat$type)) {
      warning("Skipping a pattern missing 'variable' or 'type'.")
      next
    }
    var <- pat$variable
    if (!var %in% names(data_out)) {
      warning(paste("Variable", var, "not found in data; skipping pattern."))
      next
    }
    type <- toupper(pat$type)
    if (type == "MCAR") {
      # MCAR: random missingness
      prop <- pat$prop
      if (is.null(prop) && !is.null(pat$fraction)) prop <- pat$fraction
      if (is.null(prop) && !is.null(pat$percent)) {
        prop <- pat$percent
        if (prop > 1) prop <- prop / 100
      }
      if (is.null(prop)) {
        warning(paste("MCAR pattern for", var, "has no 'prop' specified; skipping."))
        next
      }
      # Determine candidate indices (possibly restrict by region/time)
      idx_candidates <- seq_len(nrow(data_out))
      if (!is.null(pat$region) || !is.null(pat$regions)) {
        region_vals <- if (!is.null(pat$regions)) pat$regions else pat$region
        # Identify which spatial column contains these region identifiers
        spatial_cols <- names(data_out)[sapply(data_out, function(x) is.character(x) || is.factor(x))]
        for (col in spatial_cols) {
          if (any(data_out[[col]] %in% region_vals)) {
            idx_candidates <- idx_candidates[data_out[[col]][idx_candidates] %in% region_vals]
            break
          }
        }
      }
      if (!is.null(pat$start) && !is.null(pat$end) && "time" %in% names(data_out)) {
        start_date <- as.Date(pat$start)
        end_date <- as.Date(pat$end)
        idx_candidates <- idx_candidates[data_out$time[idx_candidates] >= start_date & data_out$time[idx_candidates] <= end_date]
      }
      # Randomly select a subset of these indices to set as NA
      if (length(idx_candidates) > 0) {
        n_remove <- ceiling(prop * length(idx_candidates))
        remove_idx <- sample(idx_candidates, size = n_remove, replace = FALSE)
        data_out[[var]][remove_idx] <- NA
      }
    } else if (type == "MAR") {
      # MAR: systematic missingness by region/time
      idx_sel <- seq_len(nrow(data_out))
      if (!is.null(pat$region) || !is.null(pat$regions)) {
        region_vals <- if (!is.null(pat$regions)) pat$regions else pat$region
        spatial_cols <- names(data_out)[sapply(data_out, function(x) is.character(x) || is.factor(x))]
        for (col in spatial_cols) {
          if (any(data_out[[col]] %in% region_vals)) {
            idx_sel <- idx_sel[data_out[[col]][idx_sel] %in% region_vals]
            break
          }
        }
      }
      if (!is.null(pat$start) && !is.null(pat$end) && "time" %in% names(data_out)) {
        start_date <- as.Date(pat$start)
        end_date <- as.Date(pat$end)
        idx_sel <- idx_sel[data_out$time[idx_sel] >= start_date & data_out$time[idx_sel] <= end_date]
      }
      if (length(idx_sel) > 0) {
        data_out[[var]][idx_sel] <- NA
      }
    } else if (type == "MNAR") {
      # MNAR: missingness depends on the variable's value
      if (is.null(pat$threshold) || is.null(pat$direction)) {
        warning(paste("MNAR pattern for", var, "requires 'threshold' and 'direction'; skipping."))
        next
      }
      thr <- pat$threshold
      dir <- tolower(pat$direction)
      cond_func <- NULL
      if (dir %in% c(">", "gt", "above")) {
        cond_func <- function(x) x > thr
      } else if (dir %in% c("<", "lt", "below")) {
        cond_func <- function(x) x < thr
      } else {
        warning(paste("Unrecognized direction", pat$direction, "for MNAR pattern; skipping."))
        next
      }
      idx_sel <- seq_len(nrow(data_out))
      if (!is.null(pat$region) || !is.null(pat$regions)) {
        region_vals <- if (!is.null(pat$regions)) pat$regions else pat$region
        spatial_cols <- names(data_out)[sapply(data_out, function(x) is.character(x) || is.factor(x))]
        for (col in spatial_cols) {
          if (any(data_out[[col]] %in% region_vals)) {
            idx_sel <- idx_sel[data_out[[col]][idx_sel] %in% region_vals]
            break
          }
        }
      }
      if (!is.null(pat$start) && !is.null(pat$end) && "time" %in% names(data_out)) {
        start_date <- as.Date(pat$start)
        end_date <- as.Date(pat$end)
        idx_sel <- idx_sel[data_out$time[idx_sel] >= start_date & data_out$time[idx_sel] <= end_date]
      }
      if (length(idx_sel) > 0) {
        # Apply condition to values in the subset
        vals <- data_out[[var]][idx_sel]
        to_na <- which(cond_func(vals))
        if (length(to_na) > 0) {
          data_out[[var]][idx_sel[to_na]] <- NA
        }
      }
    } else {
      warning(paste("Unknown missingness type", pat$type, "- skipping pattern."))
      next
    }
  }
  return(data_out)
}

#' Replace ideal variables with proxy variables
#' @description
#' Transform specified "ideal" predictor variables into proxy variables by applying lag, sign reversal, and/or attenuation to their values. This simulates using imperfect proxies in place of ideal predictors (for example, using a delayed or inverted indicator as a proxy). New proxy columns are added to the dataset, and (optionally) the original ideal variables can be dropped.
#'
#' @param data A data frame containing the fully observed variables (e.g., output from simulation before gaps are introduced).
#' @param proxy_map A data frame (or similar object) specifying the mapping from ideal variables to proxies and the transformations to apply. It should have at least the columns:
#' \itemize{
#'   \item \code{ideal}: Name of the ideal predictor variable in \code{data}.
#'   \item \code{proxy}: Name for the proxy variable to be created (if \code{NA} or not provided, a name will be generated by appending \code{"_proxy"} to the ideal name).
#'   \item \code{lag}: An integer number of time steps to lag the ideal variable (proxy will use the value of the ideal from \code{lag} periods earlier). Default 0 (no lag).
#'   \item \code{reverse}: Logical, if \code{TRUE}, the proxy's effect is reversed in sign (the values are multiplied by -1). Default \code{FALSE}.
#'   \item \code{attenuate}: Numeric attenuation factor for effect size. The ideal values are multiplied by this factor to produce the proxy (values between 0 and 1 will weaken the effect, values >1 would amplify it). Default 1 (no attenuation).
#' }
#' @param remove_original Logical, whether to remove the original ideal variable columns from the returned data. Default is \code{FALSE} (keep the original variables alongside their proxies).
#'
#' @details
#' For each row in \code{proxy_map}, this function creates a new proxy variable in the data. The proxy is derived from the ideal variable by:
#' (1) optionally shifting the time series by \code{lag} steps (within each region, if data contains multiple regions), so that the proxy at time \emph{t} reflects the ideal value from an earlier time \emph{t-lag};
#' (2) optionally multiplying by -1 if \code{reverse = TRUE} to invert the relationship;
#' (3) multiplying by the \code{attenuate} factor to scale the effect size.
#' Any combination of these transformations can be specified. If the data contains a \code{time} column and multiple regions, the lag is applied separately within each region's time series.
#'
#' @return A modified data frame that includes the new proxy variable columns. If \code{remove_original = TRUE}, the original ideal variables used to create proxies are dropped from the returned data.
#'
#' @examples
#' # Suppose 'health_services' is an ideal predictor in sim_data; create a lagged proxy
#' proxy_map <- data.frame(
#'   ideal = "health_services",
#'   proxy = "health_services_proxy",
#'   lag = 1,            # use previous time step's value
#'   reverse = FALSE,
#'   attenuate = 0.5     # halve the effect
#' )
#' data_with_proxy <- apply_proxies(sim_data, proxy_map, remove_original = FALSE)
#' # The new column 'health_services_proxy' is the lagged and attenuated version of 'health_services'.
#'
#' @export
apply_proxies <- function(data, proxy_map, remove_original = FALSE) {
  data_out <- data
  proxy_df <- as.data.frame(proxy_map, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(proxy_df))) {
    ideal_var <- proxy_df$ideal[i]
    if (is.null(ideal_var) || !(ideal_var %in% names(data_out))) {
      warning(paste("Ideal variable", ideal_var, "not found in data; skipping."))
      next
    }
    # Determine new proxy column name
    proxy_name <- proxy_df$proxy[i]
    if (is.null(proxy_name) || proxy_name == "" || is.na(proxy_name)) {
      proxy_name <- paste0(ideal_var, "_proxy")
    }
    if (proxy_name %in% names(data_out)) {
      warning(paste("Proxy variable name", proxy_name, "already exists in data; skipping to avoid overwrite."))
      next
    }
    # Transformation parameters
    lag_steps <- 0
    if ("lag" %in% names(proxy_df) && !is.na(proxy_df$lag[i])) {
      lag_steps <- as.integer(proxy_df$lag[i])
      if (is.na(lag_steps) || lag_steps < 0) lag_steps <- 0
    }
    reverse_flag <- FALSE
    if ("reverse" %in% names(proxy_df) && !is.na(proxy_df$reverse[i])) {
      reverse_flag <- as.logical(proxy_df$reverse[i])
    }
    attenuate_factor <- 1
    if ("attenuate" %in% names(proxy_df) && !is.na(proxy_df$attenuate[i])) {
      attenuate_factor <- as.numeric(proxy_df$attenuate[i])
    }
    # Start with ideal variable's values
    x <- data_out[[ideal_var]]
    # Apply lag if needed
    if (!is.null(lag_steps) && lag_steps > 0) {
      x_lagged <- rep(NA, length(x))
      if ("time" %in% names(data_out)) {
        # Identify lowest-level region column if present
        spatial_cols <- names(data_out)[sapply(data_out, function(z) is.character(z) || is.factor(z))]
        lowest_region_col <- NULL
        if (length(spatial_cols) > 0) {
          lowest_region_col <- spatial_cols[which.max(sapply(spatial_cols, function(c) length(unique(data_out[[c]]))))]
        }
        if (!is.null(lowest_region_col)) {
          for (reg in unique(data_out[[lowest_region_col]])) {
            idx <- which(data_out[[lowest_region_col]] == reg)
            if (length(idx) < 2) {
              # If only one observation for this region, nothing to lag
              next
            }
            if ("time" %in% names(data_out)) {
              idx <- idx[order(data_out$time[idx])]
            }
            if (lag_steps < length(idx)) {
              x_lagged[idx[(lag_steps+1):length(idx)]] <- x[idx[1:(length(idx)-lag_steps)]]
            }
            # First lag_steps entries for this region remain NA
          }
        } else {
          # No explicit region grouping column found, treat entire vector as one series
          if (lag_steps < length(x)) {
            x_lagged[(lag_steps+1):length(x)] <- x[1:(length(x)-lag_steps)]
          }
        }
      } else {
        # No time column (or single time series), apply global lag
        if (lag_steps < length(x)) {
          x_lagged[(lag_steps+1):length(x)] <- x[1:(length(x)-lag_steps)]
        }
      }
      x <- x_lagged
    }
    # Reverse effect sign if requested
    if (!is.null(reverse_flag) && reverse_flag) {
      x <- -x
    }
    # Attenuate effect size
    if (!is.null(attenuate_factor)) {
      x <- x * attenuate_factor
    }
    # Add the proxy variable to data
    data_out[[proxy_name]] <- x
  }
  if (remove_original) {
    to_drop <- unique(proxy_df$ideal)
    to_drop <- to_drop[to_drop %in% names(data_out)]
    data_out <- data_out[, !names(data_out) %in% to_drop]
  }
  return(data_out)
}

