#' Calculate reference ET using the Hargreaves method.
#'
#' @param tmin The minimum air temperature in deg C.
#' @param tmax The maximum air temperature in deg C.
#' @param tmean The average air temperature in deg C.
#' @param lat The latitude of the measurements degrees.
#' @param day The julian day of the measurement. Can either be an integer or date object.
#' @param srad Downwelling shortwave radiaiton in W m^-2. Defaults to NULL. If
#' left as null, srad will be estimated using latitude and julian day.
#'
#' @return Daily reference ET in mm day^-1
#' @export
#'
#' @examples
#' etr_hargreaves(10, 15, 5, 30, 150)
etr_hargreaves <- function(tmin, tmax, tmean, lat, day, srad = NULL) {
  checkmate::assert_multi_class(
    c(tmin, tmax, tmean, lat, day, srad),
    c("numeric", "SpatRaster")
  )
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
  } else {
    srad <- wm2_to_mj(srad)
  }
  0.0023 * (tmean + 17.8) * (tmax - tmin)**0.5 * 0.408 * srad
}
