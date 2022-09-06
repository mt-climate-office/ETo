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
  stopifnot(is.numeric(day) || inherits(day, 'Date'))
  if (inherits(day, 'Date')) {
    day <- as.numeric(format(day, '%j'))
  }
  0.409*sin((((2*pi)/365)*day)-1.39)
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
  stopifnot(is.numeric(day) || inherits(day, 'Date'))
  if (inherits(day, 'Date')) {
    day <- as.numeric(format(day, '%j'))
  }
  1+0.033*cos(((2*pi)/365)*day)
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
calc_sunset_hour_angle <- function(lat, declination) {
  acos((-tan(lat)*tan(declination)))
}

#' Calculate the theoretical extraterrestrial radiation at a given latitude.
#'
#' @param sunset_hour The sunset hour angle in radians. Can be calculated using `calc_sunset_hour_angle`
#' @param lat The latitude in radians.
#' @param declination The solar declination in radians. Can be calculated using `calc_solar_declination`
#'
#' @return
#' @export
#'
#' @examples
calc_extraterrestrial_rad <- function(sunset_hour, lat, declination) {
  ((24*60)/pi)*(0.0820*(sunset_hour*sin(lat)*sin(declination) + cos(lat)*cos(declination)*sin(sunset_hour)))
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
calc_clear_sky_radiation <- function(elev, extra_rad) {
  (0.75 + ((2 * 10^-5)* elev)) * extra_rad
}

#' Calculate the saturation vapor pressure given the daily average temperature.
#'
#' @param temp The temperature in degrees C.
#'
#' @return The saturation vapor pressure in kPa.
#' @export
#'
#' @examples
calc_sat_vapor_pressure <- function(temp) {
  0.6108*exp((17.27*temp)/(temp + 237.3))
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
calc_act_vapor_pressure <- function(svp, rh) {
  svp * (rh/100)
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
calc_radiation_fraction <- function(radiation, clear_sky) {
  frac <- radiation/clear_sky
  if (frac > 1) {
    return(1)
  } else {
    return(frac)
  }
}

#' Calculate the outgoing longwave radiation.
#'
#' @param temp The air temperature in degrees C.
#' @param act_vapor_pressure The actual vapor pressure. Can be calculated with `calc_act_vapor_pressure`.
#' @param radiation_fraction The fraction of radiation relative to clear sky. Can be calculated with `calc_radiation_fraction`.
#'
#' @return The outgoing longwave radiation in MJ m^-2 day^-1.
#' @export
#'
#' @examples
calc_longwave_radiation <- function(temp, act_vapor_pressure, radiation_fraction) {
  4.903e-9 * ((temp+273.16)^4)*(0.34-(0.14*(sqrt(act_vapor_pressure))))*((1.35*(radiation_fraction))-0.35)
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
calc_shortwave_radiation <- function(radiation, reference = 0.23) {
  (1 - 0.23) * radiation
}

#' Calculate the net outgoing radiation.
#'
#' @param shortwave The outgoing shortwave radiation in  MJ m^-2 day^-1. Can be calculated with `calc_shortwave_radiation`.
#' @param longwave The outgoing longwave radiation in x MJ m^-2 day^-1.. Can be calculated with `calc_longwave_radiation`.
#'
#' @return The net outgoing radiation in  MJ m^-2 day^-1.
#' @export
#'
#' @examples
calc_net_radiation <- function(shortwave, longwave) {
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
calc_svp_slope <- function(temp) {
  (4098 * (0.6108*exp((17.27*temp)/(temp + 237.3)))) / ((temp + 237.3)^2)
}

#' Estimate the atmospheric pressure given elevation. This assumes a constant pressure of 101.3 kPa
#'
#' @param elev The elevation in meters.
#'
#' @return The pressure in kPa.
#' @export
#'
#' @examples
calc_pressure <- function(elev) {
  # Assumes constant sea level pressure of 101.3!!!
  101.3 * (((293 - 0.0065 * elev) / 293) ^ 5.26)
}

#' Calculate the psychrometric constant given air pressure.
#'
#' @param pressure The air pressure in kPa
#'
#' @return The psychrometric constant in kPa degC^-1
#' @export
#'
#' @examples
calc_psychrometric_constant <- function(pressure) {
  0.000665 * pressure
}
