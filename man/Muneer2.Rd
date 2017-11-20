\name{Muneer2}
\alias{Muneer2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Muneer2 transposition model
}
\description{
Muneer's transposition model (Muneer et. al., 2004)
}
\usage{
Muneer2(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, DHI \eqn{D_h} and extraterrestrial GHI \eqn{I_{oh}}.}
}
\details{
Muneer2 model (Muneer et. al., 2004, Muneer, 1990) is the same as \link{Muneer1}, except that the modeling for \eqn{b} (for a sunlit surface under non-overcast sky) follows a second order polynomial whose coefficients are obtained empirically.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
T. Muneer, Solar radiation model for Europe, \emph{Building Services Engineering Research and Technology}, Volume 1, Issue 4, 1990, Pages 153-163, \url{http://dx.doi.org/10.1177/014362449001100405}. \cr
T. Muneer, C. Gueymard, H. Kambezidis, Hourly slope irradiation and illuminance, in: \emph{Solar Radiation and Daylight Models}. Butterworth-Heinemann, Oxford, 2004, Pages 143-221.
}
\seealso{
\code{\link{Muneer1}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Muneer2(NREL)

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

