# Generate Dummy Parameters for Checchi et al. (2017) DAG

Returns a data frame of plausible effect sizes and lags for use with the
Checchi et al. (2017) mortality DAG. The parameters represent estimated
contributions to crude mortality rate (CMR) in deaths per 10,000 people
per day.

## Usage

``` r
dummy_checchi_2017_parameters()
```

## Value

A data frame with columns:

- from:

  Name of the predictor node

- to:

  Name of the outcome node

- effect_size:

  Change in CMR (deaths per 10,000/day) per 1-unit increase in the
  predictor

- lag:

  Lag in days before the effect manifests

## Examples

``` r
params <- dummy_checchi_2017_parameters()
head(params)
#>                                                      from
#> 1 Exposure to armed attacks or mechanical force of nature
#> 2 Exposure to armed attacks or mechanical force of nature
#> 3 Exposure to armed attacks or mechanical force of nature
#> 4 Exposure to armed attacks or mechanical force of nature
#> 5 Exposure to armed attacks or mechanical force of nature
#> 6 Exposure to armed attacks or mechanical force of nature
#>                                           to effect_size lag
#> 1                        Forced displacement         0.2   0
#> 2          Interruption of chronic treatment         0.3  14
#> 3           Sexual and gender-based violence         0.4  14
#> 4            Burden and typology of injuries         1.2   0
#> 5 Mental health and psychosocial functioning         0.2  30
#> 6        Humanitarian public health services        -0.5   0
```
