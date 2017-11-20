\name{Muneer1}
\alias{Muneer1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Muneer1 transposition model
}
\description{
Muneer's transposition model (Muneer, 1990)
}
\usage{
Muneer1(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, DHI \eqn{D_h} and extraterrestrial GHI \eqn{I_{oh}}.}
}
\details{
Muneer's model (Muneer, 1990) treats the shaded and sunlit surfaces separately and further distinguishes between overcast and non-overcast conditions of the sunlit surface. Muneer1 uses the setting: \eqn{b = 5.73} (shaded surface); \eqn{b = 1.68} (sunlit surface under overcast sky); and \eqn{b = -0.62} (sunlit surface under non-overcast sky).

}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
T. Muneer, Solar radiation model for Europe, \emph{Building Services Engineering Research and Technology}, Volume 1, Issue 4, 1990, Pages 153-163, \url{http://dx.doi.org/10.1177/014362449001100405}.
}

\seealso{
\code{\link{Muneer2}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Muneer1(NREL)

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

