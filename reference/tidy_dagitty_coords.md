# Tidy a dagitty object with optional coordinates

Converts a dagitty DAG into a tidy data frame. Allows user-provided
coordinates, otherwise defaults to automatic layout.

## Usage

``` r
tidy_dagitty_coords(
  .dagitty,
  coords = NULL,
  seed = NULL,
  layout = "nicely",
  ...
)
```

## Arguments

- .dagitty:

  A dagitty object representing a DAG.

- coords:

  Optional. A named list of coordinates for nodes. If not provided,
  coordinates are generated.

- seed:

  Optional numeric seed for reproducible layout.

- layout:

  A layout type; defaults to `"nicely"`. `"time_ordered"` is also
  supported.

- ...:

  Additional arguments passed to the layout generation functions.

## Value

A tidy dagitty object suitable for use with ggdag functions.
