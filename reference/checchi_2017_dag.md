# Corrected Checchi et al. (2017) DAG Definition

Defines the DAG from Checchi et al. (2017), accurately including all
edges.

## Usage

``` r
checchi_2017_dag(parameters = NULL)
```

## Arguments

- parameters:

  Optional dataframe with parameters for DAG edges, with columns:

  - from: origin node

  - to: destination node

  - effect_size: numeric magnitude of effect

  - lag: numeric time lag in days (default 0)

## Value

A DAG based on Checchi et al. (2017) with class "nurah_dag".

## Examples

``` r
params <- data.frame(
  from = c("Exposure to armed attacks or mechanical force of nature", "Food insecurity",
           "Nutritional status", "Burden of endemic infectious diseases"),
  to = c("Burden and typology of injuries", "Nutritional status",
         "Burden of endemic infectious diseases", "Population mortality"),
  effect_size = c(1.0, 0.7, 1.5, 2.0),
  lag = c(0, 30, 0, 0)
)
dag <- checchi_2017_dag(parameters = params)
```
