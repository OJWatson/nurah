# Simulate a Single DAG Node

Generates simulated values for one node in the DAG, given its parent
values and predefined parameters for its distribution and relationship.

## Usage

``` r
simulate_node(node, parents_data, dag_params, noise_level = 1)
```

## Arguments

- node:

  Name of the node (used for reference in error messages or debugging).

- parents_data:

  A data frame (or list) containing the values of the parent nodes
  required to simulate this node. Each column (or list element)
  corresponds to one parent, with each row representing one observation
  (e.g., one time step for one region). If the node has no parents, this
  can be an empty data frame or list; in that case, the function will
  simulate based on baseline parameters.

- dag_params:

  A list of parameters defining how to simulate the node. Typically this
  would be the element of a DAG definition corresponding to the node. It
  should include at least:

  - `parents`: Character vector of parent node names (can be empty or
    `NULL` if none).

  - `dist`: Distribution type for the node's stochastic simulation. For
    example, `"normal"` for a continuous variable, `"poisson"` for
    counts, `"negbin"` for negative binomial, `"binary"` (or
    `"bernoulli"`) for 0/1 outcomes, or `"none"` for deterministic.

  - Additional parameters specific to the distribution and relationship.

- noise_level:

  Numeric factor to scale the stochastic noise. `noise_level = 0`
  produces deterministic output (no randomness), while `noise_level = 1`
  uses the full noise as specified in `dag_params`. Default is 1.

## Value

A numeric vector of simulated values for the node.
