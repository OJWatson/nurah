# Summarize Simulation Study Results

Aggregates results across replicates to compute mean bias, RMSE, and
coverage.

## Usage

``` r
summarize_study_results(results, by = c("scenario_id", "model_id"))
```

## Arguments

- results:

  Data frame from run_simulation_study

- by:

  Character vector of grouping variables (default: scenario_id,
  model_id)

## Value

Data frame with summary statistics by group
