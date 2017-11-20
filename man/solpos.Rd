\name{solpos}
\alias{solpos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solar positioning algorithm
}
\description{
This function performs solar positioning algorithm.
}
\usage{
solpos(Tm, lat, lon, tz, tilt, orientation, LT, alt = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tm}{Time in \code{POSIXlt} format.}
  \item{lat}{Latitude in degrees and decimal fraction.}
  \item{lon}{Longitude in degrees and decimal fraction.}
  \item{tz}{Time zone, west is negative.}
  \item{tilt}{The tilt angle of the inclined surface in degrees.}
  \item{orientation}{The azimuth angle of the inclined surface in degrees, 0 North, East positive.}
  \item{LT}{Linke turbidity. A single value (assuming the parameter is constant for a period of time), or have the same length as \code{Tm} if more granular data is available.}
  \item{alt}{Altitude in meters, default value is 0.}
}
\details{
This function calculates some frequently used solar positioning parameters. It is a wrapper for "\code{sunpos}" from the "\code{insol}" package.
}
\value{
The function returns a list containing the following solar positioning parameters.
  \item{zenith}{Zenith angle in degrees.}
  \item{incidence}{Incidence angle in degrees.}
  \item{Io}{Normal incident extraterrestrial irradiance.}
  \item{Ioh}{Extraterrestrial irradiance received by a horizontal surface.}
  \item{azimuth}{Azimuth of the Sun.}
  \item{Ics}{Ineichen-Perez clear sky GHI model.}
  \item{Icsd}{Ineichen-Perez clear sky DHI model.}
  \item{ast}{Apparent solar time. Equation of time calculation follows Lamm (1981).}
}
\references{
Ineichen, P. and Perez, R. 2002. A new airmass independent formulation for the Linke turbidity coefficient, \emph{Solar Energy}, 73(3): 151-157.\cr
Lamm, L.O. 1981. A new analytic expression for the equation of time, \emph{Solar Energy}, 26(5):465-465.
}

\section{Warning}{It is important to check whether the solar positioning is performed correctly. When the data is right averaged (averaged at the end of a period), \verb{Tm} needs to be shifted to the center of the time period.}

\examples{
library(insol)
data(NREL)
#convert datetime to POSIXct, "-30*15" is used to center the 15 min data, see Warning.
Tm <-as.POSIXct(NREL$Tm, format = "\%Y-\%m-\%d \%H:\%M:\%S", tz = "GMT")-30*15
value <- solpos(Tm, lat = 39.74, lon = -105.18, tz = -7, tilt = 40,
                orientation = 180, LT = 4.05, alt = 1829)
all.equal(value$zenith, NREL$Z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the

