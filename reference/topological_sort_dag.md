# Topologically Sort DAG Nodes

Returns a topologically sorted vector of node names from a directed
acyclic graph (DAG), ensuring that each node appears after all of its
ancestor (parent) nodes. This is useful for determining an order in
which to simulate or process nodes so that dependencies are respected.

## Usage

``` r
topological_sort_dag(dag)
```

## Arguments

- dag:

  A DAG object of class `nurah_dag` (as created by
  [`define_dag()`](https://ojwatson.github.io/nurah/reference/define_dag.md))
  containing the DAG structure and parameters.

## Value

A character vector of node names sorted in topological order (parents
always come before their children).

## Details

The function extracts the list of parent relationships from the `dag`
and then performs a topological sort using Kahn's algorithm. If the
graph is not a valid DAG (for example, if it contains a cycle or an
undefined parent reference), the function will throw an error.

All parent names must correspond to nodes defined in the DAG. If any
parent node is referenced that is not in the DAG's node set, an error is
produced.

## Examples

``` r
# Define a simple DAG: A -> B -> C
nodes <- c("A", "B", "C")
edges <- data.frame(from = c("A", "B"), to = c("B", "C"))
params <- data.frame(from = c("A", "B"), to = c("B", "C"), effect_size = c(0.5, 1.0), lag = c(0, 0))
dag_obj <- define_dag(nodes, edges, params)
topological_sort_dag(dag_obj)
#> [1] "A" "B" "C"
#> [1] "A" "B" "C"
```
