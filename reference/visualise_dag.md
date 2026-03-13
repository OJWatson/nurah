# Visualise a Directed Acyclic Graph (DAG) with Rounded‐Corner Boxes

Plots a dag, using a predefined layout if available.

## Usage

``` r
visualise_dag(
  dag,
  use_dag_layout_attr = TRUE,
  node_color = "skyblue",
  label_size = 5
)
```

## Arguments

- dag:

  A nurah_dag from
  [`define_dag()`](https://ojwatson.github.io/nurah/reference/define_dag.md)
  or
  [`checchi_2017_dag()`](https://ojwatson.github.io/nurah/reference/checchi_2017_dag.md).

- use_dag_layout_attr:

  Logical; if TRUE, use attr(dag, "layout") if present.

- node_color:

  Colour for node fill. Default \\skyblue\\.

- label_size:

  Numeric cex for labels. Default 5

## Value

Invisibly, the matrix of (x,y) coordinates used.
