\name{Temps}
\alias{Temps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Temps transposition model
}
\description{
Temps's transposition model (Temps and Coulson, 1977)
}
\usage{
Temps(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z} and incidence angle \eqn{\theta}.}
}
\details{
Temps and Coulson (1977) have suggested an anisotropic modification to the clear-sky diffuse radiance model. The additional two terms takes into account the effects of horizon brightening and circumsolar radiation.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Ralph C. Temps, K.L. Coulson, Solar radiation incident upon slopes of different orientations, \emph{Solar Energy}, Volume 19, Issue 2, 1977, Pages 179-184, ISSN 0038-092X, \url{ttp://dx.doi.org/10.1016/0038-092X(77)90056-1}.
}

\seealso{
\code{\link{Klucher}, \link{NREL}}
}
\examples{
data(NREL)
NREL$theta = NREL$theta1 #theta1 is 1st of the 5 tilted planes in the NREL dataset
NREL$s = 40 #the tilt angle corresponds to theta1
Rd <- Temps(NREL)

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
\keyword{2nd generation}
