\name{Gueymard}
\alias{Gueymard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gueymard transposition model
}
\description{
Gueymard's transposition model (Gueymard, 1987)
}
\usage{
Gueymard(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h} and DHI \eqn{D_h}.}
}
\details{
Gueymard's model (Gueymard, 1987) assumes that the radiance of a partly cloudy sky may be considered as a weighted sum of radiances of a clear and an overcast sky.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Christian Gueymard, An anisotropic solar irradiance model for tilted surfaces and its comparison with selected engineering algorithms, \emph{Solar Energy}, Volume 38, Issue 5, 1987, Pages 367-386, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(87)90009-0}. \cr
Christian Gueymard, Erratum, \emph{Solar Energy}, Volume 40, Issue 2, 1988, Page 175, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(88)90087-4}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{The original paper by Gueymard (1987) has two typographical errors. The subsequently published erratum (Gueymard, 1988) corrected the numerical values of those two coefficients. Ignoring the erratum will likely result in large prediction errors.}

\seealso{
\code{\link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Gueymard(NREL)

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

