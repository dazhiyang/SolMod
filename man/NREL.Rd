\name{NREL}
\alias{NREL}
\docType{data}
\title{
NREL irradiance data
}
\description{
A dataset containing 15 min NREL irradiance measurements and solar positioning parameters for July 2014.
}
\usage{data("NREL")}
\format{
  A data frame with 1407 observations on the following 24 variables.
  \describe{
    \item{\code{Tm}}{time (character) in "YYYY-MM-DD hh:mm:ss" format. }
    \item{\code{Gh}}{global horizontal irradiance (GHI) in \eqn{W/m^2}.}
    \item{\code{Dh}}{diffuse horizontal irradiance (GHI) in \eqn{W/m^2}.}
    \item{\code{I}}{direct normal irradiance (GHI) in \eqn{W/m^2}.}
    \item{\code{Gc1}}{global tilted irradiance (GTI) for the 40 deg South facing sensor in \eqn{W/m^2}.}
    \item{\code{Gc2}}{GTI for the 90 deg North facing sensor.}
    \item{\code{Gc3}}{GTI for the 90 deg East facing sensor.}
    \item{\code{Gc4}}{GTI for the 90 deg South facing sensor.}
    \item{\code{Gc5}}{GTI for the 90 deg West facing sensor.}
    \item{\code{rho}}{foreground's albedo (measured).}
    \item{\code{Z}}{zenith angle in degrees.}
    \item{\code{Io}}{normal incident extraterrestial irradiance in \eqn{W/m^2}.}
    \item{\code{Ioh}}{extraterrestial irradiance on a horizontal surface in \eqn{W/m^2}.}
    \item{\code{Ics}}{Ineichen-Perez clear sky GHI in \eqn{W/m^2}.}
    \item{\code{Icsd}}{Ineichen-Perez clear sky DHI in \eqn{W/m^2}.}
    \item{\code{AST}}{apparent solar time.}
    \item{\code{alpha}}{azimuth of the Sun in degrees.}
    \item{\code{theta1}}{incidence angle for the 40 deg South facing sensor in degrees.}
    \item{\code{theta2}}{incidence angle for the 90 deg North facing sensor.}
    \item{\code{theta3}}{incidence angle for the 90 deg East facing sensor.}
    \item{\code{theta4}}{incidence angle for the 90 deg South facing sensor.}
    \item{\code{theta5}}{incidence angle for the 90 deg West facing sensor.}
    \item{\code{epsilon}}{Sky's clearness.}
    \item{\code{Delta}}{Sky's brightness.}
  }
}
\details{
National Renewable Energy Laboratory (NREL), Solar Radiation Research Laboratory (SRRL), Baseline Measurement System (BMS) provides accurate 1 min weather data in Boulder, Colorado (39.74, -105.18). The site has an altitude of 1829 m, in GMT -7 time zone. The raw 1 min data is averaged into 15 min. Only filtered data from July 2014 is included in this data file.\cr

This dataset contains global tilted irradiance (GTI) data from 5 orientations. Incidence angles to these 5 surfaces, as well as other solar positioning parameters are pre-calculated using \code{\link{solpos}}.
}
\source{
\url{https://www.nrel.gov/midc/apps/go2url.pl?site=BMS}
}
\references{
A. Andreas, T. Stoffel, NREL Solar Radiation Research Laboratory (SRRL): Baseline Measurement System (BMS), Golden, Colorado (Data), NREL Report No. DA-5500-56488, 1981,
\url{http://dx.doi.org/10.5439/1052221}.
}
\examples{
data(NREL)
}
\keyword{datasets}
