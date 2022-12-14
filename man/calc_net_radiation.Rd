% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derived.R
\name{calc_net_radiation}
\alias{calc_net_radiation}
\title{Calculate the net outgoing radiation.}
\usage{
calc_net_radiation(shortwave, longwave)
}
\arguments{
\item{shortwave}{The outgoing shortwave radiation in  MJ m^-2 day^-1. Can be calculated with \code{calc_shortwave_radiation}.}

\item{longwave}{The outgoing longwave radiation in x MJ m^-2 day^-1. Can be calculated with \code{calc_longwave_radiation}.}
}
\value{
The net outgoing radiation in  MJ m^-2 day^-1.
}
\description{
Calculate the net outgoing radiation.
}
\examples{
calc_net_radiation(10, 5)
}
