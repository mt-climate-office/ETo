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
