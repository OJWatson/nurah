# Create Default Scenario Grid

Helper function to create a standard scenario grid for simulation
studies.

## Usage

``` r
create_scenario_grid(
  n_govs = 10,
  n_months = 12,
  noise_levels = c(0.5, 1),
  regimes = c("MCAR", "MAR", "MNAR"),
  pi_bases = c(0.5, 0.7, 0.9),
  mnar_omegas = c(0, -0.3, -0.6)
)
```

## Arguments

- n_govs:

  Vector of governorate counts to test

- n_months:

  Vector of month counts to test

- noise_levels:

  Vector of noise levels to test

- regimes:

  Missingness regimes to test

- pi_bases:

  Baseline observation probabilities to test

- mnar_omegas:

  MNAR sensitivity parameters to test

## Value

Data frame of scenarios with scenario_id
