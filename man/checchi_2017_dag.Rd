% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dag.R
\name{checchi_2017_dag}
\alias{checchi_2017_dag}
\title{Corrected Checchi et al. (2017) DAG Definition}
\usage{
checchi_2017_dag(parameters = NULL)
}
\arguments{
\item{parameters}{Optional dataframe with parameters for DAG edges, with columns:
\itemize{
\item from: origin node
\item to: destination node
\item effect_size: numeric magnitude of effect
\item lag: numeric time lag in days (default 0)
}}
}
\value{
A DAG based on Checchi et al. (2017) with class "nurah_dag".
}
\description{
Defines the DAG from Checchi et al. (2017), accurately including all edges.
}
\examples{
params <- data.frame(
  from = c("Exposure to armed attacks or mechanical force of nature", "Food insecurity",
           "Nutritional status", "Burden of endemic infectious diseases"),
  to = c("Burden and typology of injuries", "Nutritional status",
         "Burden of endemic infectious diseases", "Population mortality"),
  effect_size = c(1.0, 0.7, 1.5, 2.0),
  lag = c(0, 30, 0, 0)
)
dag <- checchi_2017_dag(parameters = params)
}
