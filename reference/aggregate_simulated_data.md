# Aggregate Daily Simulation Data to Periods

Aggregates detailed daily simulation outputs into the specified temporal
resolution (e.g., weekly, monthly) by region. This function sums flow
variables over each period and carries forward stock variables like
population as of the end of each period.

## Usage

``` r
aggregate_simulated_data(
  daily_data,
  time_seq_out = NULL,
  resolution = NULL,
  start_date = NULL
)
```

## Arguments

- daily_data:

  A data frame of daily data by region, containing at least a `date`
  column and various numeric indicator columns to aggregate.

- time_seq_out:

  (Optional) A Date vector of period start dates to aggregate by. If
  provided, these define the boundaries for aggregation (with an
  additional boundary at the end of the simulation). If not provided,
  `resolution` must be given.

- resolution:

  (Optional) The target time resolution for aggregation, as a string
  (e.g., `"week"`, `"month"`). Ignored if `time_seq_out` is supplied. If
  `time_seq_out` is not given, the function will derive period breaks
  from the minimum date in `daily_data` at the specified resolution.

- start_date:

  (Optional) If `time_seq_out` is not given, an anchor start date for
  the periods (defaults to the first date in `daily_data`). This date
  will be used as the start of the first period.

## Value

A data frame aggregated by region and time period. It contains one row
per region per period. The output includes all region identifier
columns, a `date` column indicating the start date of each period, and
all numeric indicator columns aggregated. Flow metrics (e.g., daily
counts) are summed over the period, while the `population` column (if
present) is taken as the last day's value in the period.
