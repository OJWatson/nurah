# Fit Model with Multiple Imputation

Uses mice package for multiple imputation, fits model to each imputed
dataset, and pools results.

## Usage

``` r
fit_with_imputation(
  data,
  final_formula,
  vars_needed,
  n_imputations = 5,
  priors = NULL,
  chains = 4,
  cores = 1,
  iter = 2000,
  warmup = floor(iter/2),
  seed = NULL,
  control = list(),
  family = gaussian(),
  ...
)
```

## Arguments

- data:

  Data frame with missing values

- final_formula:

  Model formula

- vars_needed:

  Variables needed for model

- n_imputations:

  Number of imputations

- priors:

  brms priors

- chains:

  Number of MCMC chains

- cores:

  Number of CPU cores

- iter:

  MCMC iterations

- warmup:

  Warmup iterations

- seed:

  Random seed

- control:

  Stan control parameters

- family:

  Model family

- ...:

  Additional arguments for brm

## Value

List of brmsfit objects (one per imputation)
