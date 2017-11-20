\name{Perez3}
\alias{Perez3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perez3 transposition model
}
\description{
Latest version of Perez's transposition model (Perez et al., 1990)
}
\usage{
Perez3(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, sky's brightness \eqn{\Delta} and sky's clearness \eqn{\varepsilon}.}
}
\details{
The Perez3 irradiance models offer a practical representation of solar irradiance by considering the sky hemisphere as a three-part geometrical framework, namely, the circumsolar disc (as a point source), the horizon band (as an infinitesimally thin region) and the isotropic background. Perez3 adopts the coefficients given by Perez et al. (1990).
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Richard Perez, Pierre Ineichen, Robert Seals, Joseph Michalsky, Ronald Stewart, Modeling daylight availability and irradiance components from direct and global irradiance, \emph{Solar Energy}, Volume 44, Issue 5, 1990, Pages 271-289, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(90)90055-H}.
}
\seealso{
\code{\link{Perez1}, \link{Perez2}, \link{Perez4}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Perez3(NREL)

#calculate Gc and plot
Ic <- cos(d2r(NREL$theta))/cos(d2r(NREL$Z))*(NREL$Gh-NREL$Dh)
Dg <- NREL$Gh*0.5*NREL$rho*(1-cos(d2r(NREL$s)))
Dc <- Rd*NREL$Dh
Gc <- Ic + Dg + Dc

plot(NREL$Gc1, Gc, pch = ".") #measured versus predicted
abline(0,1,col=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Transposition}
\keyword{3rd generation}

