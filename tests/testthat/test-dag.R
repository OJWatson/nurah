library(testthat)

test_that("define_dag creates correct nurah_dag object", {
  nodes <- c("Nutrition", "Food insecurity", "Mortality")
  edges <- data.frame(from = c("Food insecurity", "Nutrition"),
                      to = c("Nutrition", "Mortality"))
  parameters <- data.frame(from = c("Food insecurity", "Nutrition"),
                           to = c("Nutrition", "Mortality"),
                           effect_size = c(-0.5, 1.2),
                           lag = c(30, 0))

  dag <- define_dag(nodes, edges, parameters)

  expect_s3_class(dag, "nurah_dag")
  expect_true(inherits(dag$dag, "dagitty"))
  expect_equal(dag$parameters, parameters)
})

test_that("checchi_2017_dag returns correct DAG structure", {
  params <- data.frame(
    from = c("Exposure to armed attacks or mechanical force of nature", "Food insecurity"),
    to = c("Burden and typology of injuries", "Nutritional status"),
    effect_size = c(1.0, 0.7),
    lag = c(0, 30)
  )

  dag <- checchi_2017_dag(parameters = params)

  expect_s3_class(dag, "nurah_dag")
  expect_true(inherits(dag$dag, "dagitty"))
  expect_true(!is.null(attr(dag, "layout")))
  expect_equal(dag$parameters, params)
})

test_that("tidy_dagitty_coords correctly sets user coords", {
  nodes <- c("A", "B")
  edges <- data.frame(from = "A", to = "B")
  dagitty_obj <- define_dag(nodes, edges)$dag

  coords_df <- data.frame(name = c("A", "B"), x = c(0, 1), y = c(0, 1))
  tidy_dag <- tidy_dagitty_coords(dagitty_obj, coords = coords_df)

  expect_true("tidy_dagitty" %in% class(tidy_dag))
  expect_equal(tidy_dag$data$x, c("A" = 0, "B"= 1))
  expect_equal(tidy_dag$data$y, c("A" = 0, "B"= 1))
})

test_that("visualise_dag returns a ggplot object", {
  nodes <- c("A", "B")
  edges <- data.frame(from = "A", to = "B")
  dag <- define_dag(nodes, edges)

  plot_obj <- visualise_dag(dag)
  expect_true(inherits(plot_obj, "ggplot"))
})
