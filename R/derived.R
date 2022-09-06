calc_solar_declination <- function(day) {
  stopifnot(is.numeric(day))
  0.409*sin((((2*pi)/365)*day)-1.39)
}

calc_inverse_relative_distance <- function(day) {
  1+0.033*cos(((2*pi)/365)*day)
}

calc_sunset_hour_angle <- function(lat, dec) {
  acos((-tan(lat)*tan(dec)))
}

calc_extraterrestrial_rad <- function(sunset_hour, lat, declination) {
  ((24*60)/pi)*(0.0820*(sunset_hour*sin(lat)*sin(declination) + cos(lat)*cos(declination)*sin(sunset_hour)))
}

calc_clear_sky_radiation <- function(elev, extra_rad) {
  (0.75 + ((2 * 10^-5)* elev)) * extra_rad
}

calc_sat_vapor_pressure <- function(temp) {
  0.6108*exp((17.27*temp)/(temp + 237.3))
}

calc_act_vapor_pressure <- function(svp, rh) {
  svp * (rh/100)
}

calc_radiation_fraction <- function(radiation, clear_sky) {
  frac <- radiation/clear_sky
  if (frac > 1) {
    return(1)
  } else {
    return(frac)
  }
}

calc_longwave_radiation <- function(temp, act_vapor_pressure, radiation_fraction) {
  4.903e-9 * ((temp+273.16)^4)*(0.34-(0.14*(sqrt(act_vapor_pressure))))*((1.35*(radiation_fraction))-0.35)
}

calc_shortwave_radiation <- function(radiation, reference = 0.23) {
  (1 - 0.23) * radiation
}

calc_net_radiation <- function(shortwave, longwave) {
  shortwave - longwave
}

calc_svp_slope <- function(temp) {
  (4098 * (0.6108*exp((17.27*temp)/(temp + 237.3)))) / ((temp + 237.3)^2)
}

calc_pressure <- function(elev) {
  # Assumes constant sea level pressure of 101.3!!!
  101.3 * (((293 - 0.0065 * elev) / 293) ^ 5.26)
}

calc_psychrometric_constant <- function(pressure) {
  0.000665 * pressure
}
