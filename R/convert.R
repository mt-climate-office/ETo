#' Convert latitude in degrees to radians.
#'
#' @param lat The latitude in degrees to convert into radians.
#'
#' @return The latitude in radians.
#' @export
#'
#' @examples
#' lat_to_radians(90)
lat_to_radians <- function(lat) {
  stopifnot(is.numeric(lat), abs(lat) <= 360)
  (pi/180)*lat
}

wm2_to_mj <- function(radiation) {
  radiation * 3600 * 24 * 1e-6
}
