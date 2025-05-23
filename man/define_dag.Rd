% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dag.R
\name{define_dag}
\alias{define_dag}
\title{Define a Directed Acyclic Graph (DAG) for mortality prediction}
\usage{
define_dag(nodes, edges, parameters = NULL)
}
\arguments{
\item{nodes}{A vector of node names (health domains).}

\item{edges}{A data frame specifying directed edges between nodes with columns from and to.}

\item{parameters}{A data frame specifying parameters for each edge including:
\itemize{
\item from: node from which effect originates
\item to: node affected
\item effect_size: numeric magnitude of effect
\item lag: numeric time lag in days (default 0)
}}
}
\value{
A list containing:
\itemize{
\item dag: A dagitty object representing the DAG.
\item parameters: parameters dataframe.
}
}
\description{
Define a Directed Acyclic Graph (DAG) for mortality prediction
}
\examples{
nodes <- c("Nutrition", "Food insecurity", "Mortality")
edges <- data.frame(from = c("Food insecurity", "Nutrition"),
                    to = c("Nutrition", "Mortality"))
parameters <- data.frame(from = c("Food insecurity", "Nutrition"),
                         to = c("Nutrition", "Mortality"),
                         effect_size = c(-0.5, 1.2),
                         lag = c(30, 0))
dag <- define_dag(nodes, edges, parameters)
}
