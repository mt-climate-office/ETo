#' Get the elevation in meters for a given latitude and longitude.
#' Latitude/longitude must be in degrees/geographic projection.
#'
#' @param lat The latitude to extract elevation for. Can either be a floating point value or a
#' vector of
#' @param lon The longitude to extract elevation for.
#'
#' @return A single value
#' @export
#'
#' @examples
get_elev_from_point <- function(lat, lon) {
  data.frame(x = lon, y = lat) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 4326) |>
    elevatr::get_elev_point() %$%
    elevation
}

