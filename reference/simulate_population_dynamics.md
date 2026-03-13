# Simulate Population Dynamics with IDP Flows

Simulates changes in population over time for each region, driven by
daily inward and outward movements of internally displaced persons
(IDPs). This function updates the population for each day based on the
previous day's population and net IDP flow.

## Usage

``` r
simulate_population_dynamics(
  daily_data,
  initial_population,
  idp_in_col = "idp_in",
  idp_out_col = "idp_out"
)
```

## Arguments

- daily_data:

  A data frame of daily simulation data by region. It must contain at
  least a `date` column and one or more region identifier columns (e.g.,
  region names or codes). It should also contain columns for IDP inflows
  and outflows (or a single net flow column), unless no displacement is
  being modeled.

- initial_population:

  Numeric value or vector giving the initial population at the start of
  the simulation for each lowest-level region. If a single number is
  provided, all regions start with that population. If a vector is
  provided, it should be either:

  - A named vector, where names correspond to the region identifiers in
    `daily_data` (matching the lowest-level region names or IDs).

  - An unnamed vector with length equal to the number of lowest-level
    regions (in the same order as the regions appear when grouping
    `daily_data` by region).

- idp_in_col:

  Name of the column in `daily_data` representing IDP inflows to a
  region. Defaults to `"idp_in"`. If this column is not found, the
  function will interpret the data differently (see Details).

- idp_out_col:

  Name of the column in `daily_data` representing IDP outflows from a
  region. Defaults to `"idp_out"`. If this column is not found or set to
  `NULL`, and an inflow column is present, that inflow is assumed to
  represent net population change (positive for influx, negative for
  outflux).

## Value

The `daily_data` data frame with an additional column `population`,
giving the simulated population for each region on each day.
