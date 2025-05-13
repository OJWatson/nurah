library(testthat)

test_that("setup_spatial_structure handles numeric input correctly", {
  res <- setup_spatial_structure(3)
  expect_equal(nrow(res), 3)
  expect_named(res, c("country", "region"))
})

test_that("setup_spatial_structure handles named numeric vector", {
  res <- setup_spatial_structure(c(country=1, region=2, district=4))
  expect_equal(nrow(res), 4)
  expect_named(res, c("country", "region", "district"))
})

test_that("setup_spatial_structure handles character input", {
  res <- setup_spatial_structure(c("North", "South", "East"))
  expect_equal(nrow(res), 3)
  expect_equal(res$region, c("North", "South", "East"))
})

test_that("setup_spatial_structure returns provided dataframe as-is", {
  df <- data.frame(country="X", region=c("A","B"), district=c("A1","B1"))
  res <- setup_spatial_structure(df)
  expect_equal(res, df)
})

test_that("setup_spatial_structure errors on invalid input", {
  expect_error(setup_spatial_structure(list(a=1)), "Invalid spatial_structure specification")
})
