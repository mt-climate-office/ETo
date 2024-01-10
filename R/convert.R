#' Convert latitude in degrees to radians.
#'
#' @param lat The latitude in degrees to convert into radians.
#'
#' @return The latitude in radians.
#' @export
#'
#' @examples
#' lat_to_radians(45)
lat_to_radians <- function(lat) {
  checkmate::check_multi_class(lat, c("SpatRaster", "numeric"))
  (pi / 180) * lat
}

#' Convert radiation in watts per m^2 per day to MJ per m^2 per day.
#'
#' @param radiation Radiation in watts per m^2.
#'
#' @return Radiation in MJ per m^2.
#' @export
#'
#' @examples
#' wm2_to_mj(2)
wm2_to_mj <- function(radiation) {
  radiation * 3600 * 24 * 1e-6
}

#' Adjust wind speed measurement to a 2-meter observation.
#'
#' @param ws The wind speed in meters/second to be adjusted.
#' @param wind_height The height in meters of the wind speed observation.
#'
#' @return The wind speed in meters/second adjusted to a 2-meter observation.
#' @export
#'
#' @examples
#' adjust_wind_speed(2.777, 10)
adjust_wind_speed <- function(ws, wind_height) {
  if (wind_height != 2) {
    ws <- ws * (4.87 / log(67.8 * wind_height - 5.42))
  }

  return(ws)
}
