# Define a Directed Acyclic Graph (DAG) for mortality prediction

Define a Directed Acyclic Graph (DAG) for mortality prediction

## Usage

``` r
define_dag(nodes, edges, parameters = NULL)
```

## Arguments

- nodes:

  A vector of node names (health domains).

- edges:

  A data frame specifying directed edges between nodes with columns from
  and to.

- parameters:

  A data frame specifying parameters for each edge including:

  - from: node from which effect originates

  - to: node affected

  - effect_size: numeric magnitude of effect

  - lag: numeric time lag in days (default 0)

## Value

A list containing:

- dag: A dagitty object representing the DAG.

- parameters: parameters dataframe.

## Examples

``` r
nodes <- c("Nutrition", "Food insecurity", "Mortality")
edges <- data.frame(from = c("Food insecurity", "Nutrition"),
                    to = c("Nutrition", "Mortality"))
parameters <- data.frame(from = c("Food insecurity", "Nutrition"),
                         to = c("Nutrition", "Mortality"),
                         effect_size = c(-0.5, 1.2),
                         lag = c(30, 0))
dag <- define_dag(nodes, edges, parameters)
```
