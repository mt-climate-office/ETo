#' Calculate daily total reference ET at a given point.
#'
#' @param lat Latitude of ETo calculation in degrees.
#' @param day The julian day ETo is being calculated on.
#' @param rh Relative humidity (%).
#' @param temp Air temperature in Celsius.
#' @param srad Shortwave radiation in W m^-2.
#' @param ws Wind speed at 2 meters in M s^-1.
#' @param elev Elevation in meters.
#' @param reference The albedo of the reference surface (defaults to 0.23 for grass).
#'
#' @return Reference ET in mm per day.
#' @export
#'
#' @examples
#' etr_penman_monteith(30, 150, 88, 5, 160, 5, 1000)
etr_penman_monteith <- function(lat, day, rh, temp, srad, ws, elev, reference = 0.23) {
  checkmate::assert_multi_class(
    c(lat, day, rh, temp, srad, ws, elev, reference),
    c("numeric", "SpatRaster")
  )

  lat <- lat_to_radians(lat)
  declination <- calc_solar_declination(day)
  sunset_hour_angle <- calc_sunset_hour_angle(lat, declination)
  extraterrestrial_rad <- calc_extraterrestrial_rad(sunset_hour_angle, lat, declination)
  clear_sky_radiation <- calc_clear_sky_radiation(elev, extraterrestrial_rad)
  sat_vapor_pressure <- calc_sat_vapor_pressure(temp)
  actual_vapor_pressure <- calc_act_vapor_pressure(sat_vapor_pressure, rh)
  rad_in <- wm2_to_mj(srad)
  radiation_fraction <- calc_radiation_fraction(rad_in, clear_sky_radiation)
  rad_lw <- calc_longwave_radiation(temp, actual_vapor_pressure, radiation_fraction)
  rad_sw <- calc_shortwave_radiation(rad_in, reference)
  rad_net <- calc_net_radiation(rad_sw, rad_lw)
  svp_slope <- calc_svp_slope(temp)
  pressure <- calc_pressure(elev)
  psy <- calc_psychrometric_constant(pressure)

  ((0.408 * svp_slope * (rad_net - 0)) + (psy * (900 / (temp + 273)) * (ws * (sat_vapor_pressure - actual_vapor_pressure)))) /
    (svp_slope + psy * (1 + 0.34 * ws))
}
