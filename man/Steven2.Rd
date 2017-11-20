\name{Steven2}
\alias{Steven2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Steven2 transposition model
}
\description{
An alternative Steven's transposition model (Steven and Unsworth, 1979)
}
\usage{
Steven2(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z} and incidence angle \eqn{\theta}.}
}
\details{
Steven2 selects the model coefficients based on the closest zenith angle.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
M.D. Steven and M.H. Unsworth, The diffuse solar irradiance of slopes under cloudless skies, \emph{Quarterly Journal of the Royal Meteorological Society}, Volume 105, Issues 445, 1979, Pages 593-602, \url{http://dx.doi.org/10.1002/qj.49710544507}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Steven1}, \link{Steven3}, \link{Steven4}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Steven2(NREL)

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

