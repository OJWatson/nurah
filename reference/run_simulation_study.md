# Run Simulation Study for Mortality Estimation

Main function for running a simulation study to evaluate mortality
estimators under various scenarios (effect strengths, noise levels,
missingness regimes).

## Usage

``` r
run_simulation_study(
  dag,
  scenarios,
  models = NULL,
  n_replicates = 100,
  mortality_node = "Population mortality",
  predictors = NULL,
  mar_drivers = NULL,
  seed = NULL,
  chains = 2,
  iter = 1000,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- dag:

  A nurah_dag object defining the causal structure

- scenarios:

  Data frame defining simulation scenarios. Each row is a scenario with
  columns:

  - scenario_id: Unique identifier

  - n_govs: Number of governorates (spatial units)

  - n_months: Number of months (time periods)

  - initial_pop: Initial population per governorate

  - noise_level: Noise scaling factor (0-1)

  - missingness_regime: "MCAR", "MAR", or "MNAR"

  - pi_base: Baseline observation probability

  - mnar_omega: MNAR sensitivity parameter (if regime = "MNAR")

  - Any additional columns for scenario-specific parameters

- models:

  Named list of model specifications to fit. Each element should be a
  list with:

  - formula: Model formula (or NULL to use DAG-derived formula)

  - family: brms family (default: negbinomial())

  - complete_case: Whether to use complete case analysis (default:
    FALSE)

  - missing_method: "ipw" or "imputation" (if complete_case = FALSE)
    Default models include "complete_case", "ipw", and "full_obs"
    (oracle with no missingness)

- n_replicates:

  Number of simulation replicates per scenario (default: 100)

- mortality_node:

  Name of the mortality outcome node in the DAG

- predictors:

  Character vector of predictor variables to include in missingness

- mar_drivers:

  Character vector of variables driving MAR missingness

- seed:

  Random seed for reproducibility

- chains:

  Number of MCMC chains per model fit

- iter:

  Number of MCMC iterations per chain

- cores:

  Number of CPU cores to use

- verbose:

  Whether to print progress messages

## Value

A data frame with one row per scenario x model x replicate containing:

- scenario_id, model_id, replicate

- true_total_deaths: True total deaths from simulation

- est_total_deaths: Estimated total deaths from model

- bias: Absolute bias (estimated - true)

- rel_bias: Relative bias (bias / true)

- rmse: Root mean squared error

- coverage: Whether 95% CI contains true value

- Additional governorate-level metrics

## Examples

``` r
if (FALSE) { # \dontrun{
# Define scenarios
scenarios <- expand.grid(
  n_govs = 10,
  n_months = 12,
  initial_pop = 100000,
  noise_level = c(0.5, 1.0),
  missingness_regime = c("MCAR", "MAR", "MNAR"),
  pi_base = c(0.5, 0.8),
  mnar_omega = c(0, -0.5)
)
scenarios$scenario_id <- seq_len(nrow(scenarios))

# Run study
results <- run_simulation_study(
  dag = checchi_2017_dag(dummy_checchi_2017_parameters()),
  scenarios = scenarios,
  n_replicates = 10,
  mortality_node = "Population mortality"
)
} # }
```
