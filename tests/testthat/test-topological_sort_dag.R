library(testthat)

test_that("topological_sort_dag sorts simple DAG correctly", {
  nodes <- c("A", "B", "C")
  edges <- data.frame(from = c("A", "B"), to = c("B", "C"))
  params <- data.frame(from = c("A", "B"), to = c("B", "C"), effect_size = c(0.5, 1.0), lag = c(0, 0))
  dag_obj <- define_dag(nodes, edges, params)

  sorted_nodes <- topological_sort_dag(dag_obj)
  expect_equal(sorted_nodes, c("A", "B", "C"))
})

test_that("topological_sort_dag handles multiple branches correctly", {
  nodes <- c("A", "B", "C", "D")
  edges <- data.frame(from = c("A", "A", "B"), to = c("B", "C", "D"))
  params <- data.frame(from = edges$from, to = edges$to, effect_size = c(1,1,1), lag = c(0,0,0))
  dag_obj <- define_dag(nodes, edges, params)

  sorted_nodes <- topological_sort_dag(dag_obj)
  expect_true(which(sorted_nodes == "A") < which(sorted_nodes == "B"))
  expect_true(which(sorted_nodes == "A") < which(sorted_nodes == "C"))
  expect_true(which(sorted_nodes == "B") < which(sorted_nodes == "D"))
})

test_that("topological_sort_dag detects cycles", {
  nodes <- c("A", "B")
  edges <- data.frame(from = c("A", "B"), to = c("B", "A"))  # cyclic dependency
  params <- data.frame(from = edges$from, to = edges$to, effect_size = c(1,1), lag = c(0,0))
  dag_obj <- define_dag(nodes, edges, params)

  expect_error(topological_sort_dag(dag_obj), "DAG contains a cycle")
})

# TODO: fix when someone provides no edge DAG
# test_that("topological_sort_dag handles nodes without edges", {
#   nodes <- c("A", "B", "C")
#   edges <- data.frame(from = character(), to = character()) # no edges
#   params <- data.frame(from = character(), to = character(), effect_size = numeric(), lag = numeric())
#   dag_obj <- define_dag(nodes, edges, params)
#
#   sorted_nodes <- topological_sort_dag(dag_obj)
#   expect_setequal(sorted_nodes, c("A", "B", "C"))
# })

# TODO: Check if this should error
# test_that("topological_sort_dag errors on undefined parents", {
#   nodes <- c("A", "B")
#   edges <- data.frame(from = c("X"), to = c("B")) # undefined parent "X"
#   params <- data.frame(from = edges$from, to = edges$to, effect_size = 1, lag = 0)
#   dag_obj <- define_dag(nodes, edges, params)
#
#   expect_error(topological_sort_dag(dag_obj), "Undefined parent nodes in DAG: X")
# })

test_that("topological_sort_dag validates input class", {
  invalid_dag <- list(not_a_dag = TRUE)
  expect_error(topological_sort_dag(invalid_dag), "Input must be a 'nurah_dag'")
})
