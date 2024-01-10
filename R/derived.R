#' Given the Julian day, calculate the solar declination.
#'
#' @param day The julian day as either an integer or date.
#'
#' @return The solar declination in radians.
#' @export
#'
#' @examples
#' calc_solar_declination(60)
calc_solar_declination <- function(day) {
  checkmate::assert_multi_class(day, c("numeric", "Date", "SpatRaster"))
  if (inherits(day, "Date")) {
    day <- as.numeric(format(day, "%j"))
  }
  0.409 * sin((((2 * pi) / 365) * day) - 1.39)
}

#' Calculate the inverse relative distance between the earth and sun for a given day.
#'
#' @param day The julian day as either an integer or date.
#'
#' @return The inverse relative distance from earth to sun.
#' @export
#'
#' @examples
#' calc_inverse_relative_distance(60)
calc_inverse_relative_distance <- function(day) {
  checkmate::assert_multi_class(day, c("numeric", "Date", "SpatRaster"))

  if (inherits(day, "Date")) {
    day <- as.numeric(format(day, "%j"))
  }
  1 + 0.033 * cos(((2 * pi) / 365) * day)
}

#' Given the latitude and solar declination, calculate the sunset hour angle.
#'
#' @param lat The latitude in radians.
#' @param declination The solar declination in radians. Can be calculated using `calc_solar_declination`
#'
#' @return The solar sunset hour angle in radians
#' @export
#'
#' @examples
#' calc_sunset_hour_angle(0.5, 1)
calc_sunset_hour_angle <- function(lat, declination) {
  checkmate::assert_multi_class(lat, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(declination, c("numeric", "SpatRaster"))

  acos((-tan(lat) * tan(declination)))
}

#' Calculate the theoretical extraterrestrial radiation at a given latitude.
#'
#' @param sunset_hour The sunset hour angle in radians. Can be calculated using `calc_sunset_hour_angle`
#' @param lat The latitude in radians.
#' @param declination The solar declination in radians. Can be calculated using `calc_solar_declination`
#' @param inverse_distance The inverse relative distance between the earth and sun calculated with `calc_inverse_relative_distance`.
#'
#' @return The extraterrestrial radiation in MJ m^-2 d^-1
#' @export
#'
#' @examples
#' calc_extraterrestrial_rad(23, 0.5, 2.5, 0.96)
calc_extraterrestrial_rad <- function(sunset_hour, lat, declination, inverse_distance) {
  checkmate::assert_multi_class(sunset_hour, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(lat, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(declination, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(inverse_distance, c("numeric", "SpatRaster"))

  ((24 * 60) / pi) * (0.0820 * inverse_distance * (sunset_hour * sin(lat) * sin(declination) + cos(lat) * cos(declination) * sin(sunset_hour)))
}

#' Calculate the clear sky radiation at a location given an elevation.
#'
#' @param elev The elevation in meters.
#' @param extra_rad The extraterrestrial radiation in MJ m^-2 day^-1. Can be calculated using `calc_extraterrestrial_rad`
#'
#' @return The clear sky solar radiation in MJ m^-2 day^-1.
#' @export
#'
#' @examples
#' calc_clear_sky_radiation(1000, 10)
calc_clear_sky_radiation <- function(elev, extra_rad) {
  checkmate::assert_multi_class(elev, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(extra_rad, c("numeric", "SpatRaster"))

  (0.75 + ((2 * 10^-5) * elev)) * extra_rad
}

#' Calculate the saturation vapor pressure given the daily average temperature.
#'
#' @param temp The temperature in degrees C.
#'
#' @return The saturation vapor pressure in kPa.
#' @export
#'
#' @examples
#' calc_sat_vapor_pressure(5)
calc_sat_vapor_pressure <- function(temp) {
  checkmate::assert_multi_class(temp, c("numeric", "SpatRaster"))
  0.6108 * exp((17.27 * temp) / (temp + 237.3))
}

#' Calculate the actual vapor pressure given the saturation vapor pressure and relative humidity.
#'
#' @param svp The saturation vapor pressure in kPa. Can be calculated using `calc_sat_vapor_pressure`
#' @param rh The relative humidity (%).
#'
#' @return The actual vapor pressure in kPa.
#' @export
#'
#' @examples
#' calc_act_vapor_pressure(0.5, 88)
calc_act_vapor_pressure <- function(svp, rh) {
  checkmate::assert_multi_class(svp, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(rh, c("numeric", "SpatRaster"))
  svp * (rh / 100)
}

#' Calculate fraction of radiation relative to clear sky radiation.
#'
#' @param radiation Incoming radiation in MJ m^-2 day^-1.
#' @param clear_sky The clear sky radiation in MJ m^-2 day^-1. Can be calculated with `calc_clear_sky_radiation`.
#'
#' @return The fraction of radiation relative to clear sky (unitless).
#' @export
#'
#' @examples
#' calc_radiation_fraction(1.2, 3.5)
calc_radiation_fraction <- function(radiation, clear_sky) {
  checkmate::assert_multi_class(radiation, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(clear_sky, c("numeric", "SpatRaster"))
  frac <- radiation / clear_sky
  if (inherits(frac, "SpatRaster")) {
    frac[frac > 1] <- 1
    return(frac)
  } else {
    if (frac > 1) {
      return(1)
    } else {
      return(frac)
    }
  }
}

#' Calculate the outgoing longwave radiation.
#'
#' @param temp The air temperature in degrees C.
#' @param act_vapor_pressure The actual vapor pressure. Can be calculated with `calc_act_vapor_pressure`.
#' @param radiation_fraction The fraction of radiation relative to clear sky. Can be calculated with `calc_radiation_fraction`.
#' @param tmin The daily minimum temperature. Defaults to `NULL`. If specified with tmax, they are used in place of `temp`.
#' @param tmax The daily maximum temperature. Defaults to `NULL`. If specified with tmin, they are used in place of `temp`.
#'
#' @return The outgoing longwave radiation in MJ m^-2 day^-1.
#' @export
#'
#' @examples
#' calc_longwave_radiation(5, 1.3, 0.5)
calc_longwave_radiation <- function(temp, act_vapor_pressure, radiation_fraction, tmin = NULL, tmax = NULL) {
  checkmate::assert_multi_class(temp, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(act_vapor_pressure, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(radiation_fraction, c("numeric", "SpatRaster"))
  if (!is.null(tmin) & !is.null(tmax)) {
    temp_term <- mean(
      c(
        (tmax + 273.16)^4,
        (tmin + 273.16)^4
      )
    )
  } else {
    temp_term <- (temp + 273.16)^4
  }
  4.903e-9 * temp_term * (0.34 - (0.14 * (sqrt(act_vapor_pressure)))) * ((1.35 * (radiation_fraction)) - 0.35)
}

#' Calculate the outgoing shortwave radiation.
#'
#' @param radiation The incoming radiation in MJ m^-2 day^-1.
#' @param reference The reference surface being used for ETo calculation. (A floating point between 0 - 1)
#'
#' @return The outgoing shortwave radiation in  MJ m^-2 day^-1.
#' @export
#'
#' @examples
#' calc_shortwave_radiation(10, 0.23)
calc_shortwave_radiation <- function(radiation, reference = 0.23) {
  checkmate::assert_multi_class(radiation, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(reference, c("numeric", "SpatRaster"))

  (1 - 0.23) * radiation
}

#' Calculate the net outgoing radiation.
#'
#' @param shortwave The outgoing shortwave radiation in  MJ m^-2 day^-1. Can be calculated with `calc_shortwave_radiation`.
#' @param longwave The outgoing longwave radiation in x MJ m^-2 day^-1. Can be calculated with `calc_longwave_radiation`.
#'
#' @return The net outgoing radiation in  MJ m^-2 day^-1.
#' @export
#'
#' @examples
#' calc_net_radiation(10, 5)
calc_net_radiation <- function(shortwave, longwave) {
  checkmate::assert_multi_class(shortwave, c("numeric", "SpatRaster"))
  checkmate::assert_multi_class(longwave, c("numeric", "SpatRaster"))
  shortwave - longwave
}

#' Calculate the slope of the saturation vapor pressure curve.
#'
#' @param temp The temperature in degrees C
#'
#' @return The slope of the saturation vapor pressure curve in kPa.
#' @export
#'
#' @examples
#' calc_svp_slope(5)
calc_svp_slope <- function(temp) {
  checkmate::assert_multi_class(temp, c("numeric", "SpatRaster"))
  (4098 * (0.6108 * exp((17.27 * temp) / (temp + 237.3)))) / ((temp + 237.3)^2)
}

#' Estimate the atmospheric pressure given elevation. This assumes a constant pressure of 101.3 kPa
#'
#' @param elev The elevation in meters.
#'
#' @return The pressure in kPa.
#' @export
#'
#' @examples
#' calc_pressure(1000)
calc_pressure <- function(elev) {
  checkmate::assert_multi_class(elev, c("numeric", "SpatRaster"))
  # Assumes constant sea level pressure of 101.3!!!
  101.3 * (((293 - 0.0065 * elev) / 293)^5.26)
}

#' Calculate the psychrometric constant given air pressure.
#'
#' @param pressure The air pressure in kPa
#'
#' @return The psychrometric constant in kPa degC^-1
#' @export
#'
#' @examples
#' calc_psychrometric_constant(90)
calc_psychrometric_constant <- function(pressure) {
  checkmate::assert_multi_class(pressure, c("numeric", "SpatRaster"))
  0.000665 * pressure
}
