\name{Bugler1}
\alias{Bugler1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bugler1 transposition model
}
\description{
Bugler's transposition model (Bugler, 1977)
}
\usage{
Bugler1(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h} and DHI \eqn{D_h}.}
}
\details{
Bugler1 originates from the argument that the anisotropy in tilted diffuse irradiance should be accommodated through an additional component contributed by the circumsolar region; that additional component is 5\% of DNI: \deqn{R_d = \frac{1+\cos s}{2}+0.05\frac{I\cos \theta}{D_h}.}
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
J.W. Bugler, The determination of hourly insolation on an inclined plane using a diffuse irradiance model based on hourly measured global horizontal insolation, \emph{Solar Energy}, Volume 19, Issue 5, 1977, Pages 477-491, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(77)90103-7}.
}

\seealso{
\code{\link{Bugler2}, \link{NREL}}
}
\examples{
data(NREL)
NREL$theta = NREL$theta1 #theta1 is 1st of the 5 tilted planes in the NREL dataset
NREL$s = 40 #the tilt angle corresponds to theta1
Rd <- Bugler1(NREL)

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
