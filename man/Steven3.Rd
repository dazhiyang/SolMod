\name{Steven3}
\alias{Steven3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Steven3 transposition model
}
\description{
Steven's (isotropic) transposition model for overcast skies (Steven and Unsworth, 1980)
}
\usage{
Steven3(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z} and incidence angle \eqn{\theta}.}
}
\details{
Steven3 was developed for overcast skies. The formulation does not include the directional property. Therefore, it could be considered as a pseudo-isotropic model.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
M.D. Steven and M.H. Unsworth, The angular distribution and interception of diffuse solar radiation below overcast skies, \emph{Quarterly Journal of the Royal Meteorological Society}, Volume 106, Issues 447, 1980, Pages 57-61, \url{http://dx.doi.org/10.1002/qj.49710644705}.
}
\seealso{
\code{\link{Steven1}, \link{Steven2}, \link{Steven4}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Steven3(NREL)

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

