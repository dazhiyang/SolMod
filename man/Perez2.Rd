\name{Perez2}
\alias{Perez2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perez2 transposition model
}
\description{
Latest version of Perez's transposition model (Perez et al., 1990), with (Perez et al., 1987) coefficients.
}
\usage{
Perez2(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, DHI \eqn{D_h} and sky's brightness \eqn{\Delta}.}
}
\details{
Perez2 is same as \link{Perez3} except that Perez2 takes the set of \eqn{F_{ij}} values reported in Table 2 of Perez et al. (1987).
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Richard Perez, Pierre Ineichen, Robert Seals, Joseph Michalsky, Ronald Stewart, Modeling daylight availability and irradiance components from direct and global irradiance, \emph{Solar Energy}, Volume 44, Issue 5, 1990, Pages 271-289, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(90)90055-H}. \cr
Richard Perez, Robert Seals, Pierre Ineichen, Ronald Stewart, David Menicucci, A new simplified version of the perez diffuse irradiance model for tilted surfaces, \emph{Solar Energy}, Volume 39, Issue 3, 1987, Pages 221-231, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/S0038-092X(87)80031-2}.
}

\section{Warning}{The \eqn{\varepsilon} partitions in Perez2 are different from those in the well-known version of the Perez model, namely, \link{Perez3}. See \link{load.Perez.par} for more details.}

\seealso{
\code{\link{Perez1}, \link{Perez3}, \link{Perez4}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Perez2(NREL)

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

