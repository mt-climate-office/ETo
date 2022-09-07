#' Calculate reference ET using the Hargreaves method.
#'
#' @param tmin The minimum air temperature in deg C.
#' @param tmax The maximum air temperature in deg C.
#' @param tmean The average air temperature in deg C.
#' @param lat The latitude of the measurements in deg C.
#' @param day The julian day of the measurement. Can either be an integer or date object.
#'
#' @return Daily reference ET in mm day^-1
#' @export
#'
#' @examples
etr_hargreaves <- function(tmin, tmax, tmean, lat, day, srad = NULL) {

  if (is.null(srad)) {
    lat <- lat_to_radians(lat)
    declination <- calc_solar_declination(day)
    sunset_hour_angle <- calc_sunset_hour_angle(lat, declination)
    inverse_distance <- calc_inverse_relative_distance(day)

    # Credit to: https://github.com/woodcrafty/PyETo/blob/master/pyeto/fao.py
    tmp1 <- (24.0 * 60.0) / pi
    tmp2 <- sunset_hour_angle * sin(lat) * sin(declination)
    tmp3 <- cos(lat) * cos(declination) * sin(sunset_hour_angle)
    srad <- tmp1 * 0.0820 * inverse_distance * (tmp2 + tmp3)
  }

  0.0023 * (tmean + 17.8) * (tmax - tmin)**0.5 * 0.408 * srad
}
