\name{Willmott}
\alias{Willmott}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Willmott transposition model
}
\description{
Willmott's transposition model (Willmott, 1982)
}
\usage{
Willmott(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, DNI \eqn{I} and solar constant.}
}
\details{
Willmott's model (Willmott, 1982) is similar to \link{Hay1} and \link{Hay2}. Willmott assumes the diffuse radiation is contributed by two parts: an isotropic background and a circumsolar region. However, the isotropic diffuse radiation is approximated with a polynomial of the collector's tilt.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Cort J. Willmott, On the climatic optimization of the tilt and azimuth of flat-plate solar collectors, \emph{Solar Energy}, Volume 28, Issue 3, 1982, Pages 205-216, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(82)90159-1}.
}

\seealso{
\code{\link{Hay1}, \link{Hay2}, \link{NREL}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Willmott(NREL)

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

