\name{Hay1}
\alias{Hay1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hay1 transposition model
}
\description{
Hay's transposition model (Hay and Davies, 1980)
}
\usage{
Hay1(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, DHI \eqn{D_h} and extraterrestrial GHI \eqn{I_{oh}}.}
}
\details{
Hay's model (Hay and Davies, 1980) assumes diffuse radiation incident on a horizontal surface to be composed of circumsolar and uniform background sky-diffuse components.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
J.E. Hay and J.A. Davies, Calculation of the solar irradiance incident on an inclined surface, in: J.E. Hay and T.K. Won (Eds.), \emph{First Canadian Solar Radiation Data Workshop}, Toronto, Ontaria, Canada, Pages 59-72.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hay2}, \link{NREL}}
}

\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Hay1(NREL)

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

