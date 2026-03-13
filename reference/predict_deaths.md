# Extract Predicted Deaths from Fitted Model

Convenience function to get predicted mortality counts from a fitted
model.

## Usage

``` r
predict_deaths(fit, newdata = NULL, summary = TRUE)
```

## Arguments

- fit:

  A brmsfit object

- newdata:

  Data frame for prediction (if NULL, uses original data)

- summary:

  Whether to return summary statistics (default TRUE)

## Value

Data frame with predictions
