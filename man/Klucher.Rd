\name{Klucher}
\alias{Klucher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Klucher transposition model
}
\description{
Klucher's transposition model (Klucher, 1979)
}
\usage{
Klucher(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h} and DHI \eqn{D_h}.}
}
\details{
Klucher (1979) progressed the work of Temps and Coulson (1977) by developing an anisotropic model for all-sky conditions.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Ralph C. Temps, K.L. Coulson, Solar radiation incident upon slopes of different orientations, \emph{Solar Energy}, Volume 19, Issue 2, 1977, Pages 179-184, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(77)90056-1}. \cr
T.M. Klucher, Evaluation of models to predict insolation on tilted surfaces, \emph{Solar Energy}, Volume 23, Issue 2, 1979, Pages 111-114, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(79)90110-5}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Temps}, \link{NREL}}
}
\examples{
data(NREL)
NREL$theta = NREL$theta1 #theta1 is 1st of the 5 tilted planes in the NREL dataset
NREL$s = 40 #the tilt angle corresponds to theta1
Rd <- Klucher(NREL)

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
