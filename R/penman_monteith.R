#' Calculate daily total reference ET at a given point. Follows the steps outlined
#' in Zotarelli et al (2010) [https://unpkg.com/node-red-contrib-evapotranspiration@1.0.1/alogorithms/evapotranspiration%20ae45900.pdf](https://unpkg.com/node-red-contrib-evapotranspiration@1.0.1/alogorithms/evapotranspiration%20ae45900.pdf).
#'
#' @param lat Latitude of ETo calculation in degrees.
#' @param days The julian day ETo is being calculated on.
#' @param rh_mean Relative humidity (%).
#' @param rh_min Minumum daily relative humidity (%).
#' @param rh_max Maximum daily relative humidity (%).
#' @param t_mean Air temperature in Celsius.
#' @param t_min Minimum daily air temperature in degC.
#' @param t_max Maximum daily air temperature in degC.
#' @param srad Shortwave radiation in W m^-2.
#' @param ws Wind speed at 2 meters in M s^-1.
#' @param elev Elevation in meters.
#' @param reference The albedo of the reference surface (defaults to 0.23 for grass).
#' @param wind_height The height of the wind observation in meters. (defaults to 2). If it is not 2 meters, the wind speed will be corrected to a 2 meter observation.
#'
#' @return Reference ET in mm per day.
#' @export
#'
#' @examples
#' etr_penman_monteith(
#'   t_mean = NULL,
#'   t_max = 21.5,
#'   t_min = 12.3,
#'   rh_mean = NULL,
#'   rh_max = 84,
#'   rh_min = 63,
#'   lat = 50.8,
#'   days = 187,
#'   ws = 2.77778,
#'   wind_height = 10,
#'   elev = 100,
#'   reference = 0.23,
#'   srad = 255.4398
#' )
etr_penman_monteith <- function(
    lat, days, srad, ws, elev, rh_mean = NULL, rh_min = NULL, rh_max = NULL, t_mean = NULL, t_min = NULL, t_max = NULL, reference = 0.23, wind_height = 2) {
  # checkmate::assert_multi_class(
  #   c(lat, days, rh_mean, t_mean, srad, ws, elev, reference),
  #   c("numeric", "SpatRaster")
  # )

  if (!is.null(t_mean)) {
    if (!is.null(t_min) | !is.null(t_max)) {
      stop("If t_mean is not NA, neither t_min nor t_max may be passed as arguments")
    }
  }

  if (!is.null(rh_mean)) {
    if (!is.null(rh_min) | !is.null(rh_max)) {
      stop("If rh_mean is not NA, neither rh_min nor rh_max may be passed as arguments")
    }
  }

  if (!is.null(t_min) | !is.null(t_max)) {
    if (is.null(t_min) | is.null(t_max)) {
      stop("If either t_min or t_max is passed as an argument, the other must also be used.")
    }
    if (!is.null(t_mean)) {
      warning("Because t_min and t_max are both included, t_mean will be ignored.")
    }

    # Step 1: Calculate Mean Daily Temperature
    t_mean <- (t_min + t_max)/2
  }

  if (!is.null(rh_min) | !is.null(rh_max)) {
    if (is.null(rh_min) | is.null(rh_max)) {
      stop("If either rh_min or rh_max is passed as an argument, the other must also be used.")
    }
    if (!is.null(rh_mean)) {
      warning("Because rh_min and rh_max are both included, rh_mean will be ignored.")
    }
    rh_mean <- (rh_min + rh_max)/2
  }

  # Step 2: Mean Daily Solar Radiation
  # Should be 22.07 in FAO daily ETo example.
  rad_in <- wm2_to_mj(srad) # 22.07

  # Step 3: Wind Speed.
  # Adjust wind speed if observation is not at 2m above surface.
  # Should correct to 2.078 in the FAO ETo example.
  if (wind_height != 2) {
    ws <- ws * (4.87 / log(67.8 * wind_height - 5.42))
  }

  # Step 4: Slope of saturation vapor pressure curve.
  # Should be 0.122 in the FAO ETo example.
  svp_slope <- calc_svp_slope(t_mean)

  # Step 5: Calculate atmospheric pressure.
  # Should be 100.1 in the FAO ETo example.
  pressure <- calc_pressure(elev)

  # Step 6: Calculate psychrometric constant
  # Should be 0.0666 in the FAO ETo example.
  psy <- calc_psychrometric_constant(pressure)

  # Steps 7, 8 and 9 are not included here, as they are just subsets of the
  # final ETo equation.

  # Step 10: Calculate saturation vapor pressure.
  # Should be 1.997 with FAO example.
  if (!is.null(t_min) & !is.null(t_max)) {
    sat_vp_tmin <- calc_sat_vapor_pressure(t_min)
    sat_vp_tmax <- calc_sat_vapor_pressure(t_max)
    sat_vapor_pressure <- (sat_vp_tmin + sat_vp_tmax) / 2
  } else {
    sat_vapor_pressure <- calc_sat_vapor_pressure(t_mean)
  }

  # Step 11: Calculate actual vapor pressure.
  # Should be 1.409 in FAO example.
  if (!is.null(rh_min) & !is.null(rh_max)) {
    # Use all min/max values.
    if (!is.null(t_min) & !is.null(t_max)) {
      actual_vapor_pressure <- (
        calc_act_vapor_pressure(sat_vp_tmin, rh_max) +
        calc_act_vapor_pressure(sat_vp_tmax, rh_min)
      ) / 2
      # Use daily mean rh_mean and temperature.
    } else {
      actual_vapor_pressure <- calc_act_vapor_pressure(sat_vapor_pressure, rh_mean)
    }
    # If you don't have daily min/max rh_mean.
  } else {
    if (!is.null(t_min) & !is.null(t_max)) {
      actual_vapor_pressure <- calc_act_vapor_pressure(
        (sat_vp_tmin + sat_vp_tmax) / 2, rh_mean
      )
    } else {
      actual_vapor_pressure <- calc_act_vapor_pressure(sat_vapor_pressure, rh_mean)
    }
  }

  # Step 12: Calculate inverse relative Earth-Sun distance and solar declination.
  declination <- calc_solar_declination(days)
  inverse_distance <- calc_inverse_relative_distance(days)

  # Step 13: Convert latitude from degrees to radians.
  lat <- lat_to_radians(lat)

  # Step 14: Sunset hour angle.
  sunset_hour_angle <- calc_sunset_hour_angle(lat, declination)
  # Step 15: Extraterrestrial radiation
  # Should be 41.09 in the FAO ETo example.
  extraterrestrial_rad <- calc_extraterrestrial_rad(
    sunset_hour_angle, lat, declination, inverse_distance
  )

  # Step 16: Clear sky solar radiation
  # Should be 30.90 in the FAO ETo example.
  clear_sky_radiation <- calc_clear_sky_radiation(elev, extraterrestrial_rad) # 30.90

  # Step 17: Net shortwave radiation
  # Should be 17.0 in the FAO ETo example.
  rad_sw <- calc_shortwave_radiation(rad_in, reference)

  # Step 18: Net outgoing long wave radiation
  # Should be 0.71 in the FAO ETo example.
  radiation_fraction <- calc_radiation_fraction(rad_in, clear_sky_radiation)
  # Should be 3.71 in the FAO ETo example.
  rad_lw <- calc_longwave_radiation(t_mean, actual_vapor_pressure, radiation_fraction, t_min, t_max)

  # Step 19: Calculate net radiation
  rad_net <- calc_net_radiation(rad_sw, rad_lw)

  # Final Step: Calculate ETo
  eto = ((0.408 * svp_slope * (rad_net - 0)) + (psy * (900 / (t_mean + 273)) * (ws * (sat_vapor_pressure - actual_vapor_pressure)))) /
    (svp_slope + psy * (1 + 0.34 * ws))

  eto[eto < 0] <- 0
  return(eto)
}
