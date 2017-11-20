\name{Liu}
\alias{Liu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Liu transposition model
}
\description{
Liu and Jordan's isotropic transposition model (Liu and Jordan, 1961)
}
\usage{
Liu(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}}
}
\details{
An isotropic transposition model assumes a constant radiance sky. Liu and Jordan is the earliest and simplest transposition model; it is often used as a reference model.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
B.Y.H. Liu, R.C. Jordan, Daily insolation on surfaces tilted towards the equator. \emph{ASHRAE Journal}, Volume 3, 1961, Pages 53-59.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Badescu}, \link{Koronakis}, \link{Tian}}
}
\examples{
data <- data.frame(s = 40)
Rd <- Liu(data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Transposition}
\keyword{1st generation}
