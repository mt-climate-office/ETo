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
#' xs = c(-113.994, -111.032)
#' ys = c(46.8721, 45.6815)
#'
#' # Get elevation of Missoula and Bozeman.
#' elev <- get_elev_from_point(ys, xs)
#' }
get_elev_from_point <- function(lat, lon) {
  elev <- data.frame(x = lon, y = lat) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 4326) |>
    elevatr::get_elev_point()

  elev$elevation
}

#' Download a raster of elevation that matches the shape, resolution and
#' projection of an input raster.
#'
#' @param r A `terra::rast` specifying the domain to download data for.
#' @param z The zoom level download elevation data at ranging from 1 - 14.
#' One if coarser resolution, 14 is finer resolution. For more information,
#' look at `?elevatr::get_elev_raster`
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' r <- terra::rast(rh)[[1]]
#' elev <- get_elev_from_raster(r, 9)
#' terra::plot(elev)
#' }
get_elev_from_raster <- function(r, z) {
  checkmate::assert_class(r, "SpatRaster")
  checkmate::assert_number(z)

  terra::boundaries(r) |>
    terra::as.polygons() |>
    sf::st_as_sf() |>
    elevatr::get_elev_raster(z = z, clip="bbox") |>
    terra::rast() |>
    terra::project(r)
}

#' Calculate spatial daily reference ET using either the Penman-Monteith
#' method or the Hargreaves method.
#'
#' @param tmean A `terra::rast` of daily average temperature in degrees C.
#' @param tmin A `terra::rast` of  daily min temperature in degrees C.
#' @param tmax A `terra::rast` of  daily max temperature in degrees C.
#' @param srad A `terra::rast` of  daily downwelling shortwave solar radiaiton in W m^-2
#' @param rh A `terra::rast` of  daily relative humidity (%).
#' @param ws A `terra::rast` of  daily average wind speed at 2m in M s^-1
#' @param elev A `terra::rast` of elevation in meters across the domain. If the
#' argument is left blank, elevation data are downloaded using the `elevatr`
#' package.
#' @param day The Julian day of the year ETo is being calculated for. Can either
#' be an integer or a date object.
#' @param lat A `terra::rast` of latitude with the same resolution, projection,
#' and extent as all other input rasters. If left blank, a latitude raster will
#' be derived from the `tmean` input.
#' @param reference The albedo of the reference surface ranging from 0 - 1.
#' Defaults to 0.23 for grass.
#' @param z The zoom level download elevation data at ranging from 1 - 14.
#' One if coarser resolution, 14 is finer resolution. For more information,
#' look at `?elevatr::get_elev_raster`
#' @param method The method to calculate ETo. Can either be "penman" or "hargreaves".
#'
#' @return A `terra::rast` of daily reference ET.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Calculate ETo across Montana for 2015-01-01
#' srad <- terra::rast(srad) %>% terra::subset(1)
#' tmean <- terra::rast(tmean) %>% terra::subset(1)
#' tmax <- terra::rast(tmax) %>% terra::subset(1)
#' tmin <- terra::rast(tmin) %>% terra::subset(1)
#' rh <- terra::rast(rh) %>% terra::subset(1)
#' ws <- terra::rast(ws) %>% terra::subset(1)
#'
#' penman <- calc_etr_spatial(
#'   tmean = tmean, srad = srad, rh = rh, ws = ws,
#'   method = "penman", reference = 0.23, z = 9
#' )
#' hargreaves <- calc_etr_spatial(
#'   tmean = tmean, tmax = tmax, tmin = tmin, method = "hargreaves"
#' )
#' }
calc_etr_spatial <- function(
  tmean, tmin = NULL, tmax = NULL, srad = NULL, rh = NULL, ws = NULL, elev = NULL,
  day = NULL, lat = NULL, reference = 0.23, z = 9, method = "penman"
) {

  checkmate::assert_choice(method, c("penman", "hargreaves"))
  if (is.null(elev)) {
    elev <- get_elev_from_raster(tmean, z = z)
  }

  if (is.null(day)) {
    day <- terra::time(tmean)
  }

  if (is.null(lat)) {
    lat <- terra::deepcopy(tmean)
    lat[] <- terra::xyFromCell(lat, 1:terra::ncell(lat))[,2]
  }

  if (method == "penman") {
    checkmate::assert_class(tmean, "SpatRaster")
    checkmate::assert_class(srad, "SpatRaster")
    checkmate::assert_class(rh, "SpatRaster")
    checkmate::assert_class(ws, "SpatRaster")

    eto <- etr_penman_monteith(
      lat = lat, day = day, rh = rh, temp = tmean, rad = srad, ws = ws,
      elev = elev, reference = reference
    )
  } else {
    checkmate::assert_class(tmean, "SpatRaster")
    checkmate::assert_class(tmax, "SpatRaster")
    checkmate::assert_class(tmin, "SpatRaster")

    eto <- etr_hargreaves(
      tmin = tmin, tmax = tmax, tmean = tmean, lat = lat, day = day
    )
  }

  return(eto)

}

#' Calculate a timeseries of ETo across a set of input rasters. Assumes that all
#' inputs share the same resolution, extent, projection and number of layers.
#'
#' @param tmean A multilayer `terra::rast`
#' @param tmin
#' @param tmax
#' @param srad
#' @param rh
#' @param ws
#' @param elev
#' @param days
#' @param reference
#' @param z
#' @param method
#'
#' @return
#' @export
#'
#' @examples
calc_etr_spatial_stacked <- function(
    tmean, tmin = NULL, tmax = NULL, srad = NULL, rh = NULL, ws = NULL, elev = NULL,
    days = NULL, reference = 0.23, z = 9, method = "penman"
) {

  # tmean %<>% terra::rast()
  # tmean = tmean - 273.15
  # tmin %<>% terra::rast()
  # tmin = tmin - 273.15
  # tmax %<>% terra::rast()
  # tmax = tmax - 273.15
  # srad %<>% terra::rast()
  # rh %<>% terra::rast()
  # ws %<>% terra::rast()
  # elev = NULL
  # days = NULL
  # reference = 0.23
  # z = 9

  if (is.null(elev)) {
    elev <- get_elev_from_raster(tmean[[1]], z = z)
  }

  if (is.null(days)) {
    days <- terra::time(tmean) |>
      lapply(function(x) {
        jday <- format(x, "%j") |>
          as.numeric()
        temp <- terra::deepcopy(tmean[[1]])
        temp[] <- jday
        terra::time(temp) <- x
        temp
      }) |>
      terra::rast()
  }

  lat <- terra::deepcopy(tmean[[1]])
  lat[] <- terra::xyFromCell(lat, 1:terra::ncell(lat))[,2]

  ref <- terra::deepcopy(tmean[[1]])
  ref[] <- reference

  if (method == "penman") {
    checkmate::assert_class(tmean, "SpatRaster")
    checkmate::assert_class(srad, "SpatRaster")
    checkmate::assert_class(rh, "SpatRaster")
    checkmate::assert_class(ws, "SpatRaster")


    dataset <- terra::sds(lat, days, rh, tmean, srad, ws, elev, ref)

    v_etr <- Vectorize(etr_penman_monteith)
    ETo <- terra::lapp(x = dataset, v_etr)

  } else {

    dataset <- terra::sds(tmin, tmax, tmean, lat, days) %>%
      magrittr::set_names(c("tmin", "tmax", "tmean", "lat", "days", "srad"))

    v_hg <- Vectorize(etr_hargreaves)
    ETo <- terra::lapp(dataset, v_hg)

  }
}
# list.files("~/MCO_onedrive/General/nexgddp_cmip6_montana/data-derived/nexgddp_cmip6/",
#            full.names = T, pattern = "MRI-ESM2-0") %>%
#   grep("ssp585", ., value = T) %>%
#   grep(".json", ., value = T, invert = T) -> f_list
#
#
# rh <- f_list[1] %>% terra::rast() %>% terra::subset(1:365) %>% terra::wrap()
# srad <- f_list[5] %>% terra::rast() %>% terra::subset(1:365)%>% terra::wrap()
# ws <- f_list[6] %>% terra::rast()%>% terra::subset(1:365) %>% terra::wrap()
# tmean <- f_list[7] %>% terra::rast()%>% terra::subset(1:365) %>% terra::wrap()
# tmax <- f_list[8] %>% terra::rast()%>% terra::subset(1:365)%>% terra::wrap()
# tmin <- f_list[9] %>% terra::rast() %>% terra::subset(1:365)%>% terra::wrap()
# use_data(rh, srad, ws, tmean, tmax, tmin, overwrite = T)
