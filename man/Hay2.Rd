\name{Hay2}
\alias{Hay2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hay2 transposition model
}
\description{
An alternative version of Hay's transposition model (Hay, 1993)
}
\usage{
Hay2(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, DNI \eqn{I} and solar constant.}
}
\details{
This is an alternative version of \link{Hay1}. The "anisotropic index" (\eqn{A_I}) is calculated using DNI and solar constant.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
J.E. Hay, Calculating solar radiation for inclined surfaces: Pratical approaches, \emph{Renewable Energy}, Volume 3, Issues 4-5, 1993, Pages 373-380, ISSN 0960-1481, \url{http://dx.doi.org/10.1016/0960-1481(93)90104-O}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hay1}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Hay2(NREL)

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

