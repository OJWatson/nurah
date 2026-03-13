# Simulate Mortality Using Negative Binomial Distribution

Simulates mortality counts Y_it ~ NegBin(mu_it, phi) with: log(mu_it) =
log(N_it) + eta_it eta_it = intercept + sum_j beta_j \* parent_j(it) +
random effects (optional)

## Usage

``` r
simulate_mortality_negbin(
  daily_data,
  dag,
  mortality_node,
  phi = 1,
  intercept = -8,
  noise_level = 1,
  n_days,
  n_regions
)
```

## Arguments

- daily_data:

  Daily simulation data with population column

- dag:

  DAG object containing parameters

- mortality_node:

  Name of mortality outcome node

- phi:

  Overdispersion parameter (higher = less overdispersion; phi -\> Inf
  gives Poisson)

- intercept:

  Baseline log-rate intercept

- noise_level:

  Noise scaling factor

- n_days:

  Number of days in simulation

- n_regions:

  Number of regions

## Value

daily_data with mortality column added
