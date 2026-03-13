# Compute Recovery Metrics for Simulation Study

Computes bias, RMSE, and coverage metrics comparing estimated to true
values.

## Usage

``` r
compute_recovery_metrics(
  true_value,
  est_value,
  est_sd = NULL,
  ci_lower = NULL,
  ci_upper = NULL
)
```

## Arguments

- true_value:

  True value from simulation

- est_value:

  Estimated value from model

- est_sd:

  Standard deviation of estimate

- ci_lower:

  Lower bound of confidence interval

- ci_upper:

  Upper bound of confidence interval

## Value

List with bias, rel_bias, rmse, coverage
