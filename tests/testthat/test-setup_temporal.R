test_that("setup_temporal_structure handles end_date correctly", {
  res <- setup_temporal_structure("2025-01-01", end_date="2025-01-10", resolution="day")
  expect_equal(length(res$time_seq_daily), 10)
  expect_equal(res$time_seq_daily[1], as.Date("2025-01-01"))
  expect_equal(res$time_seq_daily[10], as.Date("2025-01-10"))
})

test_that("setup_temporal_structure handles n_periods correctly", {
  res <- setup_temporal_structure("2025-01-01", n_periods=2, resolution="month")
  expect_equal(res$time_seq_out, as.Date(c("2025-01-01", "2025-02-01")))
  expect_equal(res$time_seq_daily[1], as.Date("2025-01-01"))
  expect_equal(res$time_seq_daily[length(res$time_seq_daily)], as.Date("2025-02-28"))
})

test_that("setup_temporal_structure gives error if both or neither end_date/n_periods", {
  expect_error(setup_temporal_structure("2025-01-01", resolution="day"),
               "Please provide either end_date or n_periods")
  expect_error(setup_temporal_structure("2025-01-01", "2025-02-01", 3, "month"),
               "Please provide either end_date or n_periods")
})

test_that("setup_temporal_structure handles non-period-aligned end_date", {
  res <- setup_temporal_structure("2025-01-01", end_date="2025-02-10", resolution="month")
  expect_true(as.Date("2025-02-01") %in% res$time_seq_out)
  expect_true(as.Date("2025-03-01") %in% res$time_seq_out) # Should cover partial month
})

test_that("setup_temporal_structure returns correct types", {
  res <- setup_temporal_structure("2025-01-01", n_periods=1, resolution="week")
  expect_type(res$time_seq_daily, "double")
  expect_type(res$time_seq_out, "double")
  expect_s3_class(res$time_seq_daily, "Date")
  expect_s3_class(res$time_seq_out, "Date")
})
