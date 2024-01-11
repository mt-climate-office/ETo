#' Get the elevation in meters for a given latitude and longitude.
#' Latitude/longitude must be in degrees with a geographic projection.
#'
#' @param lat The latitude to extract elevation for. Can either be a floating point value or a
#' vector of latitudes.
#' @param lon The longitude to extract elevation for. Can either be a floating point value or a
#' vector of latitudes.
#'
#' @return A numeric giving the elevation in meters for the input lat(s) and lon(s).
#' @export
#'
#' @examples
#' \dontrun{
#' xs <- c(-113.994, -111.032)
#' ys <- c(46.8721, 45.6815)
#'
#' # Get elevation of Missoula and Bozeman.
#' elev <- get_elev_from_point(ys, xs)
#' }
get_elev_from_point <- function(lat, lon) {
  elev <- data.frame(x = lon, y = lat) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 4326) |>
    elevatr::get_elev_point()

  return(elev$elevation)
}

#' Download a raster of elevation that matches the shape, resolution and
#' projection of an input raster.
#'
#' @param r A `terra::rast` specifying the domain to download data for.
#' @param z The zoom level download elevation data at ranging from 1 - 14.
#' One if coarser resolution, 14 is finer resolution. For more information,
#' look at `?elevatr::get_elev_raster`
#' @param verbose Whether to display all messages when getting elevation data.
#' Defaults to FALSE.
#'
#' @return `terra::rast` of elevation in m for the input domain
#' @export
#'
#' @examples
#' \dontrun{
#' # test
#' r <- terra::rast(rh)[[1]]
#' elev <- get_elev_from_raster(r, 9)
#' terra::plot(elev)
#' }
get_elev_from_raster <- function(r, z, verbose = FALSE) {
  checkmate::assert_class(r, "SpatRaster")
  checkmate::assert_number(z)

  r |>
    terra::rast() |>
    terra::as.polygons() |>
    sf::st_as_sf() |>
    elevatr::get_elev_raster(z = z, verbose = verbose) |>
    terra::rast() |>
    terra::project(r)
}

#' Create a `terra::rast()` of julian days from a timeseries input.
#'
#' @param r The timeseries `terra::rast()` used as an input.
#'
#' @return A `terra::rast()` that is of size `terra::nlyr(r)`, where each layer is a static value representing the julian day.
#' @export
#'
#' @examples
#' \dontrun{
#' r <- terra::rast(rh)
#' get_days_from_raster(r)
#' }
get_days_from_raster <- function(r) {
  terra::time(r) |>
    lapply(function(x) {
      jday <- format(x, "%j") |>
        as.numeric()
      temp <- terra::deepcopy(r[[1]])
      temp[] <- jday
      terra::time(temp) <- x
      temp
    }) |>
    terra::rast()
}

#' Create a latitude raster from an input raster template.
#'
#' @param r The `terra::rast` to generate a latitude grid for.
#'
#' @return A single layer `terra::rast` where each pixel's value is the latitude.
#' @export
#'
#' @examples
#' @examples
#' \dontrun{
#' r <- terra::rast(rh)
#' get_lat_from_raster(r)
#' }
get_lat_from_raster <- function(r) {
  lat <- terra::deepcopy(r[[1]])
  lat[] <- terra::xyFromCell(lat, 1:terra::ncell(lat))[, 2]
  return(lat)
}

#' Calculate a timeseries of ETo across a set of input rasters.
#'
#' @description  Assumes that all inputs share the same resolution, extent, projection. Also assumes that
#' rasters have the number of layers with each layer corresponding to a day.
#' For hargreaves, only t_mean, t_min and t_max are required arguments. For PM,
#' t_mean, srad, rh, and ws must all be supplied.
#'
#' @param t_mean A multilayer timeseries `terra::rast` of mean daily temperature in degrees C.
#' @param t_min A multilayer timeseries `terra::rast` of max daily temperature in degrees C.
#' @param t_max A multilayer timeseries `terra::rast` of min daily temperature in degrees C.
#' @param srad A multilayer timeseries `terra::rast` of daily solar radiation in W m^-2.
#' @param rh A multilayer timeseries `terra::rast` of mean daily relative humidity (%).
#' @param rh_min A multilayer timeseries `terra::rast` of minimum daily relative humidity (%; defaults to NULL).
#' @param rh_max A multilayer timeseries `terra::rast` of maximum daily relative humidity (%).
#' @param ws A multilayer timeseries `terra::rast` of mean daily wind speed at 2m height in m s^-1.
#' @param elev A `terra::rast` of elevation in meters. If left blank, elevation will
#' be derived using the `elevatr` package.
#' @param days A multilayer timeseries `terra::rast` where each day is a constant value
#' of the Julian day. If left blank, it is assumed the `t_mean` input has a time
#' attribute and the raster will be derived using `terra::time`.
#' @param reference The albedo of the reference surface ranging from 0 - 1.
#' Defaults to 0.23 for grass.
#' @param z The zoom level download elevation data at ranging from 1 - 14.
#' One if coarser resolution, 14 is finer resolution. For more information,
#' look at `?elevatr::get_elev_raster`
#' @param wind_height The height of the wind observation in meters. (defaults to 2).
#' If it is not 2 meters, the wind speed will be corrected to a 2 meter observation.
#' @param method The method to calculate ETo. Can either be "penman" or "hargreaves".
#'
#' @return A `terra::rast` timeseries of ETo for the input domain.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load data. Need to read with terra::rast to unpack to a raster.
#' srad <- terra::rast(srad) |> terra::subset(1:10)
#' t_mean <- terra::rast(t_mean) |> terra::subset(1:10)
#' # Convert from K to C
#' t_mean <- t_mean - 273.15
#' t_max <- terra::rast(t_max) |> terra::subset(1:10)
#' # Convert from K to C
#' t_max <- t_max - 273.15
#' t_min <- terra::rast(t_min) |> terra::subset(1:10)
#' # Convert from K to C
#' t_min <- t_min - 273.15
#' rh <- terra::rast(rh) |> terra::subset(1:10)
#' ws <- terra::rast(ws) |> terra::subset(1:10)
#'
#' penman <- calc_etr_spatial(
#'   t_mean = t_mean, srad = srad, rh = rh, ws = ws,
#'   method = "penman", reference = 0.23, z = 3
#' )
#' hargreaves <- calc_etr_spatial(
#'   t_mean = t_mean, t_max = t_max, t_min = t_min, method = "hargreaves", z = 3
#' )
#' }
calc_etr_spatial <-
  function(
      t_mean = NULL,
      t_min = NULL,
      t_max = NULL,
      srad = NULL,
      rh = NULL,
      rh_min = NULL,
      rh_max = NULL,
      ws = NULL,
      elev = NULL,
      days = NULL,
      reference = 0.23,
      z = 9,
      wind_height = 2,
      method = "penman") {
    checkmate::assert_choice(method, c("penman", "hargreaves"))

    if(!is.null(srad)){
      model <- srad
    }else{
      model <- t_mean
    }
    if (is.null(elev)) {
        elev <- get_elev_from_raster(model[[1]], z = z)
    }

    if (is.null(days)) {
      days <- get_days_from_raster(model)
    }

    lat <- get_lat_from_raster(model)

    if (method == "penman") {

      ref <- terra::deepcopy(model[[1]])
      ref[] <- reference

      checkmate::assert_multi_class(t_mean, c("SpatRaster", "NULL"))
      checkmate::assert_multi_class(t_min, c("SpatRaster", "NULL"))
      checkmate::assert_multi_class(t_max, c("SpatRaster", "NULL"))
      checkmate::assert_class(srad, "SpatRaster")
      checkmate::assert_multi_class(rh, c("SpatRaster", "NULL"))
      checkmate::assert_multi_class(rh_min, c("SpatRaster", "NULL"))
      checkmate::assert_multi_class(rh_max, c("SpatRaster", "NULL"))
      checkmate::assert_class(ws, "SpatRaster")

      # dataset <- terra::sds(lat, days, rh, t_mean, srad, ws, elev, ref)

      ETo <- etr_penman_monteith(
        lat = lat,
        days = days,
        rh_mean = rh,
        rh_min = rh_min,
        rh_max = rh_max,
        t_mean = t_mean,
        t_min = t_min,
        t_max = t_max,
        srad = srad,
        ws = ws,
        elev = elev,
        reference = reference,
        wind_height = wind_height
      )
      #     # v_etr <- Vectorize(etr_penman_monteith)
      #     ETo <- terra::lapp(x = dataset, etr_penman_monteith, usenames = TRUE)
    } else {
      # dataset <- terra::sds(t_min, t_max, t_mean, lat, days)
      # v_hg <- Vectorize(etr_hargreaves)
      # ETo <- terra::lapp(dataset, v_hg)
      ETo <- etr_hargreaves(tmin = t_min, tmax = t_max, tmean = t_mean, lat = lat, days = days)
    }

    terra::time(ETo) <- terra::time(days)
    return(ETo)
  }
