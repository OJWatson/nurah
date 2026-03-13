# Fit DAG-based Bayesian Hierarchical Model

This function constructs and fits a Bayesian hierarchical linear model
(via **brms** and Stan) to estimate "indirect mortality" given a
user-specified causal DAG and observed data. It uses the DAG to
determine which predictors (direct causes of the outcome) to include,
and adds specified spatial hierarchy levels as random intercept effects.

## Usage

``` r
fit_dag(
  data,
  dag,
  outcome,
  spatial_levels = NULL,
  complete_case = FALSE,
  missing_method = c("ipw", "imputation"),
  ipw_formula = NULL,
  ipw_stabilize = TRUE,
  n_imputations = 5,
  priors = NULL,
  chains = 4,
  cores = getOption("mc.cores", 1),
  iter = 2000,
  warmup = floor(iter/2),
  seed = NULL,
  control = list(),
  formula_override = NULL,
  family = gaussian(),
  offset_col = NULL,
  ...
)
```

## Arguments

- data:

  A data frame of observed data. Must contain the `outcome` variable and
  all predictor and grouping variables.

- dag:

  A `nurah_dag` object defining the causal structure. This DAG is used
  to extract the direct causes of `outcome` as fixed-effect predictors.

- outcome:

  Name (string) of the dependent variable in `data` to model.

- spatial_levels:

  Character vector of spatial grouping variables for hierarchical random
  effects.

- complete_case:

  Logical. If TRUE, drop rows with any missing values (default behavior
  pre-update). If FALSE (default), handle missingness using IPW or
  imputation.

- missing_method:

  Method for handling missing data when complete_case = FALSE. Options:
  "ipw" (inverse probability weighting), "imputation" (multiple
  imputation via mice). Default is "ipw".

- ipw_formula:

  Formula for the observation/missingness model (logistic regression
  predicting whether each row is fully observed). If NULL and
  missing_method = "ipw", a default formula using available non-outcome
  columns will be constructed.

- ipw_stabilize:

  Logical. If TRUE, use stabilized IPW weights (recommended). Default
  TRUE.

- n_imputations:

  Number of imputations if missing_method = "imputation". Default 5.

- priors:

  Optional brms prior specification. If `NULL`, defaults to weakly
  informative priors.

- chains:

  Number of MCMC chains to run (defaults to 4).

- cores:

  Number of CPU cores to use for parallel sampling.

- iter:

  Number of iterations per chain (including warmup; default 2000).

- warmup:

  Number of warmup (burn-in) iterations per chain.

- seed:

  Random seed for reproducibility.

- control:

  A list of control parameters passed to Stan.

- formula_override:

  Optional formula object to use as the model formula.

- family:

  The likelihood family for the outcome. Defaults to
  [`gaussian()`](https://rdrr.io/r/stats/family.html) for continuous
  outcomes. For mortality counts, use `negbinomial()` with offset.

- offset_col:

  Name of the population/offset column for count models. If provided,
  log(offset) will be added to the model.

- ...:

  Additional arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

A `brmsfit` object (or list of fits for imputation) containing the
fitted Bayesian model.

## Details

**Missingness Handling (NEW):** When `complete_case = FALSE`, the
function implements one of two approaches:

1.  **IPW (Inverse Probability Weighting)**: Estimates P(fully observed
    \| covariates) using logistic regression, then weights each
    observation by 1/P in the outcome model. This corrects for selection
    bias under MAR (missing at random) assumptions. Stabilized weights
    (recommended) multiply by marginal P(observed) to reduce variance.

2.  **Multiple Imputation**: Uses the `mice` package to generate
    multiple imputed datasets, fits the model to each, and pools results
    using Rubin's rules.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit with IPW to handle missingness
fit <- fit_dag(
  data = sim_data_missing,
  dag = dag,
  outcome = "Population mortality",
  spatial_levels = "region",
  complete_case = FALSE,
  missing_method = "ipw",
  family = negbinomial(),
  offset_col = "population"
)
} # }
```
