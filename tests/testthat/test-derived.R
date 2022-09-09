test_that("Solar declination works as expected", {
  expect_equal(calc_solar_declination(1), -0.4010081, tolerance = 1e-3)
})

test_that("Bad type for solar declination throws error.", {
  expect_error(calc_solar_declination('a'))
})

test_that("Good types work as expected for inverse relative distance", {
  expect_equal(
    calc_inverse_relative_distance(as.Date("2022-01-01")), 1.032995, tolerance = 1e-3
  )
  expect_equal(
    calc_inverse_relative_distance(1), 1.032995, tolerance = 1e-3
  )
})

test_that("Bad type for inverse relative distance throws error.", {
  expect_error(calc_inverse_relative_distance("2022-01-01"))
})

test_that("Sunset hour works as expected", {
  expect_equal(calc_sunset_hour_angle(pi, -0.4), 1.570796, tolerance = 1e-3)
})

test_that("Bad type for sunset hour throws error.", {
  expect_error(calc_sunset_hour_angle('a', 1))
  expect_error(calc_sunset_hour_angle(1, 'a'))
})

test_that("Extraterrestrial radiation works as expected", {
  expect_equal(calc_extraterrestrial_rad(0.2, pi, -0.4), -6.8777, tolerance = 1e-3)
})

test_that("Bad type for extraterrestrial radiation throws error.", {
  expect_error(calc_extraterrestrial_rad(0.2, pi, 'a'))
  expect_error(calc_extraterrestrial_rad('a', pi, -0.4))
  expect_error(calc_extraterrestrial_rad(0.2, 'a', -0.4))
})

test_that("Clear sky radiation works as expected", {
  expect_equal(calc_clear_sky_radiation(1000, 10), 7.7, tolerance = 1e-3)
})

test_that("Bad type for clear sky radiation throws error.", {
  expect_error(calc_clear_sky_radiation(1000, 'a'))
  expect_error(calc_clear_sky_radiation('a', 1000))
})

test_that("Sat vapor pressure works as expected", {
  expect_equal(calc_sat_vapor_pressure(10),1.227963, tolerance = 1e-3)
})

test_that("Bad type for sat vapor pressure throws error.", {
  expect_error(calc_sat_vapor_pressure('a'))
})

test_that("Act vapor pressure works as expected", {
  expect_equal(calc_act_vapor_pressure(1.22, 88), 1.0736, tolerance = 1e-3)
})

test_that("Bad type for act vapor pressure throws error.", {
  expect_error(calc_act_vapor_pressure('a', 88))
  expect_error(calc_act_vapor_pressure(1.22, 'a'))
})

test_that("Radiation fraction works as expected", {
  expect_equal(calc_radiation_fraction(1.22, 88), 0.01386364, tolerance = 1e-3)
})

test_that("Bad type for radiation fraction throws error.", {
  expect_error(calc_radiation_fraction(1.22, 'a'))
  expect_error(calc_radiation_fraction('a', 88))
})

test_that("Longwave radiation works as expected", {
  expect_equal(calc_longwave_radiation(5, 1.22, 0.5), 1.768282, tolerance = 1e-3)
})

test_that("Bad type for longwave radiation throws error.", {
  expect_error(calc_longwave_radiation('a', 1.22, 0.5))
  expect_error(calc_longwave_radiation(5, 'a', 0.5))
  expect_error(calc_longwave_radiation(5, 1.22, 'a'))
})

test_that("Shortwave radiation works as expected", {
  expect_equal(calc_shortwave_radiation(1.22, 0.23), 0.9394, tolerance = 1e-3)
})

test_that("Bad type for shortwave radiation throws error.", {
  expect_error(calc_shortwave_radiation(1.22, 'a'))
  expect_error(calc_shortwave_radiation('a', 0.23))
})

test_that("Net radiation works as expected", {
  expect_equal(calc_net_radiation(1, 2), -1, tolerance = 1e-3)
})

test_that("Bad type for net radiation throws error.", {
  expect_error(calc_net_radiation(1.22, 'a'))
  expect_error(calc_net_radiation('a', 0.23))
})

test_that("SVP slope works as expected", {
  expect_equal(calc_svp_slope(5), 0.06088867, tolerance = 1e-3)
})

test_that("Bad type for svp slope throws error.", {
  expect_error(calc_svp_slope('a'))
})

test_that("Pressure works as expected", {
  expect_equal(calc_pressure(1000), 90.02462, tolerance = 1e-3)
})

test_that("Bad type for pressure throws error.", {
  expect_error(calc_pressure('a'))
})

test_that("Psychrometric constant works as expected", {
  expect_equal(calc_psychrometric_constant(90), 0.05985, tolerance = 1e-3)
})

test_that("Bad type for psychrometric constant throws error.", {
  expect_error(calc_psychrometric_constant('a'))
})

