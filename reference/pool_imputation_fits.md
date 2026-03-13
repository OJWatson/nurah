# Pool Results from Multiple Imputation Fits

Combines posterior samples from multiple brmsfit objects using Rubin's
rules.

## Usage

``` r
pool_imputation_fits(fits)
```

## Arguments

- fits:

  List of brmsfit objects from fit_with_imputation

## Value

Data frame with pooled coefficient estimates, SE, and confidence
intervals
