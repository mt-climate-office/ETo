test_that("Penman monteith works as expected.", {
  expect_equal(
    etr_penman_monteith(30, 150, 50, 5, 15, 2, 100, 0.25),
    1.539795, tolerance = 1e-4
  )
})
