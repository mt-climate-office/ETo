test_that("Getting point elevation works.", {
  elev <- get_elev_from_point(47, -115)
  expect_equal(elev, 1459.06)
})

test_that("Getting raster elevation works.", {
  elev <- get_elev_from_raster(terra::rast(rh), 2)
  checkmate::expect_class(elev, "SpatRaster")
})


test_that("Spatial multi-band ETo calculation works.", {
  tmean %<>% terra::rast() %>% terra::subset(1:10)
  tmean <- tmean - 273.15
  tmin %<>% terra::rast() %>% terra::subset(1:10)
  tmin <- tmin - 273.15
  tmax %<>% terra::rast() %>% terra::subset(1:10)
  tmax <- tmax - 273.15
  srad %<>% terra::rast() %>% terra::subset(1:10)
  rh %<>% terra::rast() %>% terra::subset(1:10)
  ws %<>% terra::rast() %>% terra::subset(1:10)

  pm <- calc_etr_spatial(
    tmean = tmean, srad = srad, rh = rh, ws = ws, z = 3, method = "penman"
  )

  har <- calc_etr_spatial(
    tmean = tmean, tmin = tmin, tmax = tmax, z = 3, method = "hargreaves"
  )

  checkmate::expect_class(pm, "SpatRaster")
  checkmate::expect_class(har, "SpatRaster")

  expect_equal(terra::nlyr(pm), 10)
  expect_equal(terra::nlyr(har), 10)
})
