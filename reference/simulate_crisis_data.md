# Simulate Crisis Data Based on a DAG Model

Generates a synthetic crisis dataset over time and across regions, using
a directed acyclic graph (DAG) of relationships between indicators. This
high-level function orchestrates the simulation by setting up time and
space dimensions, simulating each indicator (node) in the DAG, handling
dynamic population changes due to IDP flows, and aggregating results to
the desired time resolution.

## Usage

``` r
simulate_crisis_data(
  dag = NULL,
  start_date,
  end_date = NULL,
  n_periods = NULL,
  resolution = "month",
  spatial_structure = 1,
  initial_population = 10000,
  noise_level = 1,
  mortality_node = NULL,
  mortality_phi = 1,
  mortality_intercept = -8
)
```

## Arguments

- dag:

  A DAG object of class `nurah_dag` defining the nodes (indicators) and
  their relationships. Typically created by
  [`define_dag()`](https://ojwatson.github.io/nurah/reference/define_dag.md)
  or a specific DAG constructor (e.g.,
  [`checchi_2017_dag()`](https://ojwatson.github.io/nurah/reference/checchi_2017_dag.md)).
  If `NULL`, a default DAG will be used (if available).

- start_date:

  Start date of the simulation (as a `Date` or string in "YYYY-MM-DD"
  format).

- end_date:

  End date of the simulation (as a `Date` or string). Either `end_date`
  or `n_periods` must be provided (but not both).

- n_periods:

  Number of time periods to simulate at the specified output resolution
  (ignored if `end_date` is provided).

- resolution:

  Time resolution for the output data. One of `"day"`, `"week"`,
  `"month"`, `"quarter"`, or `"year"`. Defaults to `"month"`.

- spatial_structure:

  Specification of the spatial layout (regions and optional hierarchy).
  Can be:

  - A data frame with one row per lowest-level region (and optional
    higher level columns for hierarchy).

  - A named numeric vector defining counts at each hierarchy level (top
    to bottom). For example, `c(country = 1, region = 2, district = 4)`
    means 1 country, 2 regions total, 4 districts total (distributed as
    evenly as possible under the regions).

  - A single numeric value (number of lowest-level regions; one
    top-level grouping will be added if \>1 region is implied).

  - A character vector of region names (defining the lowest level; one
    top-level grouping will be added since no higher levels are
    specified).

  Defaults to `1` (a single region with one top-level grouping if
  needed).

- initial_population:

  Initial population at the start date for each lowest-level region. Can
  be a single number (applied to all regions) or a vector. If a vector,
  it can be:

  - A named vector, where names correspond to region identifiers (at the
    lowest level).

  - An unnamed vector with length equal to the number of lowest-level
    regions (in the same order as the regions in the spatial structure).

  Defaults to 10000 for all regions.

- noise_level:

  Numeric factor to scale stochastic noise in node simulations.
  `noise_level = 0` produces fully deterministic output, while
  `noise_level = 1` uses the full specified randomness. Defaults to 1.

- mortality_node:

  Character name of the mortality outcome node. If specified, this node
  will be simulated using negative binomial distribution with population
  offset.

- mortality_phi:

  Overdispersion parameter for negative binomial mortality. Default 1.

- mortality_intercept:

  Baseline log-rate for mortality (on log scale). Default -8
  (approximately 0.03 deaths per 1000 per day).

## Value

A data frame of simulated indicators aggregated at the requested time
resolution (default: monthly). Contains one row per governorate per
month with columns for each simulated indicator including mortality
counts (Y_it) and population (N_it).

## Examples

``` r
# Simulate 10 governorates x 12 months with mortality outcome
dag <- checchi_2017_dag(parameters = dummy_checchi_2017_parameters())
sim_data <- simulate_crisis_data(
  dag,
  start_date = "2023-10-01",
  n_periods = 12,
  resolution = "month",
  spatial_structure = paste0("gov", 1:10),
  initial_population = 100000,
  noise_level = 1,
  mortality_node = "Population mortality"
)
head(sim_data)
#>    country region       date gov_id time_id
#> 1 country1   gov1 2023-10-01     31      31
#> 2 country1   gov1 2023-11-01     30      60
#> 3 country1   gov1 2023-12-01     31      93
#> 4 country1   gov1 2024-01-01     31     124
#> 5 country1   gov1 2024-02-01     29     145
#> 6 country1   gov1 2024-03-01     31     186
#>   Exposure to armed attacks or mechanical force of nature
#> 1                                                       0
#> 2                                                       0
#> 3                                                       0
#> 4                                                       0
#> 5                                                       0
#> 6                                                       0
#>   Burden and typology of injuries Forced displacement Food insecurity
#> 1                               0                   0               0
#> 2                               0                   0               0
#> 3                               0                   0               0
#> 4                               0                   0               0
#> 5                               0                   0               0
#> 6                               0                   0               0
#>   Feeding and care practices Humanitarian public health services
#> 1                          0                                   0
#> 2                          0                                   0
#> 3                          0                                   0
#> 4                          0                                   0
#> 5                          0                                   0
#> 6                          0                                   0
#>   Interruption of chronic treatment Nutritional status
#> 1                                 0                  0
#> 2                                 0                  0
#> 3                                 0                  0
#> 4                                 0                  0
#> 5                                 0                  0
#> 6                                 0                  0
#>   Burden of endemic infectious diseases Epidemic occurrence and severity
#> 1                                     0                                0
#> 2                                     0                                0
#> 3                                     0                                0
#> 4                                     0                                0
#> 5                                     0                                0
#> 6                                     0                                0
#>   Reproductive and neonatal health Sexual and gender-based violence
#> 1                                0                                0
#> 2                                0                                0
#> 3                                0                                0
#> 4                                0                                0
#> 5                                0                                0
#> 6                                0                                0
#>   Mental health and psychosocial functioning Addiction Burden of NCDs
#> 1                                          0         0              0
#> 2                                          0         0              0
#> 3                                          0         0              0
#> 4                                          0         0              0
#> 5                                          0         0              0
#> 6                                          0         0              0
#>   Population mortality population
#> 1                  731      1e+05
#> 2                  830      1e+05
#> 3                  893      1e+05
#> 4                  768      1e+05
#> 5                  836      1e+05
#> 6                  564      1e+05
```
