\name{Steven4}
\alias{Steven4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Steven4 transposition model
}
\description{
A newly proposed Steven's transposition model which considers sky conditions based on sky's clearness \eqn{\varepsilon}.
}
\usage{
Steven4(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta} and sky's clearness \eqn{\varepsilon}.}
}
\details{
This model uses \link{Steven2} for \eqn{\epsilon < 1.1}, \link{Steven3} otherwise.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
M.D. Steven and M.H. Unsworth, The diffuse solar irradiance of slopes under cloudless skies, \emph{Quarterly Journal of the Royal Meteorological Society}, Volume 105, Issues 445, 1979, Pages 593-602, \url{http://dx.doi.org/10.1002/qj.49710544507}. \cr
M.D. Steven and M.H. Unsworth, The angular distribution and interception of diffuse solar radiation below overcast skies, \emph{Quarterly Journal of the Royal Meteorological Society}, Volume 106, Issues 447, 1980, Pages 57-61, \url{http://dx.doi.org/10.1002/qj.49710644705}.
}

\seealso{
\code{\link{Steven1}, \link{Steven2}, \link{Steven3}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Steven4(NREL)

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

