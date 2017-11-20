\name{Perez4}
\alias{Perez4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perez4 transposition model
}
\description{
1987 version of Perez's transposition model (Perez et al., 1987), with (Perez et al., 1990) coefficients.
}
\usage{
Perez4(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, sky's brightness \eqn{\Delta} and sky's clearness \eqn{\varepsilon}.}
}
\details{
Perez4 is same as \link{Perez1} except that the model coefficients and \eqn{\varepsilon} partitions from (Perez et al., 1990) are used.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Richard Perez, Pierre Ineichen, Robert Seals, Joseph Michalsky, Ronald Stewart, Modeling daylight availability and irradiance components from direct and global irradiance, \emph{Solar Energy}, Volume 44, Issue 5, 1990, Pages 271-289, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(90)90055-H}. \cr
Richard Perez, Robert Seals, Pierre Ineichen, Ronald Stewart, David Menicucci, A new simplified version of the perez diffuse irradiance model for tilted surfaces, \emph{Solar Energy}, Volume 39, Issue 3, 1987, Pages 221-231, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/S0038-092X(87)80031-2}.
}

\seealso{
\code{\link{Perez1}, \link{Perez2}, \link{Perez3}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Perez4(NREL)

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

