\name{Perez1}
\alias{Perez1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perez1 transposition model
}
\description{
1987 version of Perez's transposition model (Perez et al., 1987).
}
\usage{
Perez1(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, DHI \eqn{D_h} and sky's brightness \eqn{\Delta}.}
}
\details{
Perez1 assumes that the energy from the horizon band is coming from an infinitesimally thin region, but the modeling of the circumsolar region follows the original Perez model (Perez et al., 1986). Coefficients in Table 1 of (Perez et al., 1987) are used.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Richard Perez, Robert Seals, Pierre Ineichen, Ronald Stewart, David Menicucci, A new simplified version of the perez diffuse irradiance model for tilted surfaces, \emph{Solar Energy}, Volume 39, Issue 3, 1987, Pages 221-231, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/S0038-092X(87)80031-2}.\cr
R. Perez, R. Stewart, C. Arbogast, R. Seals, J. Scott, An anisotropic hourly diffuse radiation model for sloping surfaces: Description, performance validation, site dependency evaluation, \emph{Solar Energy}, Volume 36, Issue 6, 1986, Pages 481-497, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(86)90013-7}.
}

\section{Warning}{The \eqn{\varepsilon} partitions in Perez1 are different from those in the well-known version of the Perez model, namely, \link{Perez3}. See \link{load.Perez.par} for more details.}

\seealso{
\code{\link{Perez2}, \link{Perez3}, \link{Perez4}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Perez1(NREL)

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

