% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dag.R
\name{topological_sort_dag}
\alias{topological_sort_dag}
\title{Topologically Sort DAG Nodes}
\usage{
topological_sort_dag(dag)
}
\arguments{
\item{dag}{A DAG object of class \code{nurah_dag} (as created by \code{define_dag()})
containing the DAG structure and parameters.}
}
\value{
A character vector of node names sorted in topological order (parents always
come before their children).
}
\description{
Returns a topologically sorted vector of node names from a directed acyclic graph
(DAG), ensuring that each node appears after all of its ancestor (parent) nodes.
This is useful for determining an order in which to simulate or process nodes so
that dependencies are respected.
}
\details{
The function extracts the list of parent relationships from the \code{dag}
and then performs a topological sort using Kahn's algorithm. If the graph is not a
valid DAG (for example, if it contains a cycle or an undefined parent reference),
the function will throw an error.

All parent names must correspond to nodes defined in the DAG. If any parent node is
referenced that is not in the DAG's node set, an error is produced.
}
\examples{
# Define a simple DAG: A -> B -> C
nodes <- c("A", "B", "C")
edges <- data.frame(from = c("A", "B"), to = c("B", "C"))
params <- data.frame(from = c("A", "B"), to = c("B", "C"), effect_size = c(0.5, 1.0), lag = c(0, 0))
dag_obj <- define_dag(nodes, edges, params)
topological_sort_dag(dag_obj)
#> [1] "A" "B" "C"
}
