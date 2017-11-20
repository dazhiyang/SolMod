\name{Skartveit}
\alias{Skartveit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Skartveit transposition model
}
\description{
Skartveit and Olseth's transposition model (Skartveit and Olseth, 1986)
}
\usage{
Skartveit(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, DHI \eqn{D_h} and extraterrestrial GHI \eqn{I_{oh}}.}
}
\details{
Skartveit and Olseth's transposition model (Skartveit and Olseth, 1986) is modified on top of \link{Hay1}. Their model expresses diffuse irradiance as the sum of the circumsolar irradiance, the zenith sky irradiance, and the isotropic sky irradiance at unobscured horizon.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Arvid Skartveit, Jan Asle Olseth, A model for the diffuse fraction of hourly global radiation, \emph{Solar Energy}, Volume 38, Issue 4, 1987, Pages 271-274, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(87)90049-1}.
}


\seealso{
\code{\link{Hay1}, \link{Hay2}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Skartveit(NREL)

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

