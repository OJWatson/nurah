# Compute Inverse Probability Weights for Missingness

Fits a logistic regression model for P(fully observed \| covariates) and
returns inverse probability weights.

## Usage

``` r
compute_ipw_weights(
  data,
  vars_needed,
  outcome,
  ipw_formula = NULL,
  stabilize = TRUE
)
```

## Arguments

- data:

  Data frame

- vars_needed:

  Variables needed for the outcome model

- outcome:

  Name of outcome variable

- ipw_formula:

  Optional formula for observation model

- stabilize:

  Whether to use stabilized weights

## Value

List with data (complete cases only) and weights vector
