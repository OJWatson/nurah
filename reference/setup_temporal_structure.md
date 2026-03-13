# Set up Temporal Simulation Structure

Determines the sequence of daily dates and the sequence of output period
start dates based on a start date and either an end date or a number of
periods, plus a temporal resolution.

## Usage

``` r
setup_temporal_structure(
  start_date,
  end_date = NULL,
  n_periods = NULL,
  resolution = c("day", "week", "month", "quarter", "year")
)
```

## Arguments

- start_date:

  Date (or string in "YYYY-MM-DD" format) for the start of the
  simulation.

- end_date:

  Date (or string) for the end of the simulation. Either `end_date` or
  `n_periods` must be provided (but not both).

- n_periods:

  Integer number of time periods to simulate at the output resolution
  (ignored if `end_date` is provided).

- resolution:

  Time resolution for output periods. Options include `"day"`, `"week"`,
  `"month"`, `"quarter"`, or `"year"`.

## Value

A list with two Date vectors: `time_seq_daily` (all daily dates in the
simulation) and `time_seq_out` (the start date of each output period at
the specified resolution).

## Details

If `end_date` is provided, the daily sequence runs from `start_date` to
`end_date` inclusive. If `n_periods` is provided instead, the function
generates that many consecutive periods (at the given resolution)
starting from `start_date`, and computes an `end_date` to cover all
those periods fully. The output vector `time_seq_out` contains the start
date of each period.

For example, if `start_date = "2025-01-15"` and `resolution = "month"`
with `n_periods = 2`, the output periods will start on 2025-01-15 and
2025-02-15, and the daily sequence will run until 2025-03-14 to cover
two full monthly periods (since the second period runs from 2025-02-15
to 2025-03-14).

## Examples

``` r
# Using an end date
setup_temporal_structure(start_date = "2025-01-01", end_date = "2025-01-10", resolution = "day")
#> $time_seq_daily
#>  [1] "2025-01-01" "2025-01-02" "2025-01-03" "2025-01-04" "2025-01-05"
#>  [6] "2025-01-06" "2025-01-07" "2025-01-08" "2025-01-09" "2025-01-10"
#> 
#> $time_seq_out
#>  [1] "2025-01-01" "2025-01-02" "2025-01-03" "2025-01-04" "2025-01-05"
#>  [6] "2025-01-06" "2025-01-07" "2025-01-08" "2025-01-09" "2025-01-10"
#> 
setup_temporal_structure(start_date = "2025-01-01", end_date = "2025-03-10", resolution = "month")
#> $time_seq_daily
#>  [1] "2025-01-01" "2025-01-02" "2025-01-03" "2025-01-04" "2025-01-05"
#>  [6] "2025-01-06" "2025-01-07" "2025-01-08" "2025-01-09" "2025-01-10"
#> [11] "2025-01-11" "2025-01-12" "2025-01-13" "2025-01-14" "2025-01-15"
#> [16] "2025-01-16" "2025-01-17" "2025-01-18" "2025-01-19" "2025-01-20"
#> [21] "2025-01-21" "2025-01-22" "2025-01-23" "2025-01-24" "2025-01-25"
#> [26] "2025-01-26" "2025-01-27" "2025-01-28" "2025-01-29" "2025-01-30"
#> [31] "2025-01-31" "2025-02-01" "2025-02-02" "2025-02-03" "2025-02-04"
#> [36] "2025-02-05" "2025-02-06" "2025-02-07" "2025-02-08" "2025-02-09"
#> [41] "2025-02-10" "2025-02-11" "2025-02-12" "2025-02-13" "2025-02-14"
#> [46] "2025-02-15" "2025-02-16" "2025-02-17" "2025-02-18" "2025-02-19"
#> [51] "2025-02-20" "2025-02-21" "2025-02-22" "2025-02-23" "2025-02-24"
#> [56] "2025-02-25" "2025-02-26" "2025-02-27" "2025-02-28" "2025-03-01"
#> [61] "2025-03-02" "2025-03-03" "2025-03-04" "2025-03-05" "2025-03-06"
#> [66] "2025-03-07" "2025-03-08" "2025-03-09" "2025-03-10"
#> 
#> $time_seq_out
#> [1] "2025-01-01" "2025-02-01" "2025-03-01" "2025-04-01"
#> 

# Using number of periods (e.g., 3 monthly periods starting Jan 1, 2025)
setup_temporal_structure(start_date = "2025-01-01", n_periods = 3, resolution = "month")
#> $time_seq_daily
#>  [1] "2025-01-01" "2025-01-02" "2025-01-03" "2025-01-04" "2025-01-05"
#>  [6] "2025-01-06" "2025-01-07" "2025-01-08" "2025-01-09" "2025-01-10"
#> [11] "2025-01-11" "2025-01-12" "2025-01-13" "2025-01-14" "2025-01-15"
#> [16] "2025-01-16" "2025-01-17" "2025-01-18" "2025-01-19" "2025-01-20"
#> [21] "2025-01-21" "2025-01-22" "2025-01-23" "2025-01-24" "2025-01-25"
#> [26] "2025-01-26" "2025-01-27" "2025-01-28" "2025-01-29" "2025-01-30"
#> [31] "2025-01-31" "2025-02-01" "2025-02-02" "2025-02-03" "2025-02-04"
#> [36] "2025-02-05" "2025-02-06" "2025-02-07" "2025-02-08" "2025-02-09"
#> [41] "2025-02-10" "2025-02-11" "2025-02-12" "2025-02-13" "2025-02-14"
#> [46] "2025-02-15" "2025-02-16" "2025-02-17" "2025-02-18" "2025-02-19"
#> [51] "2025-02-20" "2025-02-21" "2025-02-22" "2025-02-23" "2025-02-24"
#> [56] "2025-02-25" "2025-02-26" "2025-02-27" "2025-02-28" "2025-03-01"
#> [61] "2025-03-02" "2025-03-03" "2025-03-04" "2025-03-05" "2025-03-06"
#> [66] "2025-03-07" "2025-03-08" "2025-03-09" "2025-03-10" "2025-03-11"
#> [71] "2025-03-12" "2025-03-13" "2025-03-14" "2025-03-15" "2025-03-16"
#> [76] "2025-03-17" "2025-03-18" "2025-03-19" "2025-03-20" "2025-03-21"
#> [81] "2025-03-22" "2025-03-23" "2025-03-24" "2025-03-25" "2025-03-26"
#> [86] "2025-03-27" "2025-03-28" "2025-03-29" "2025-03-30" "2025-03-31"
#> 
#> $time_seq_out
#> [1] "2025-01-01" "2025-02-01" "2025-03-01"
#> 
```
