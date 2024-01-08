test_that("Penman monteith works as expected.", {
  expect_equal(
    etr_penman_monteith(
      t_mean = NULL,
      t_max = 21.5,
      t_min = 12.3,
      rh_mean = NULL,
      rh_max = 84,
      rh_min = 63,
      lat = 50.8,
      days = 187,
      ws = 2.77778,
      wind_height = 10,
      elev = 100,
      reference = 0.23,
      srad = 255.4398
    ),
    3.88004,
    tolerance = 1e-4
  )
})
