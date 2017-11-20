\name{Badescu}
\alias{Badescu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Badescu transposition model
}
\description{
Badescu's transposition model (Badescu, 2002)
}
\usage{
Badescu(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}}
}
\details{
This is a pseudo-isotropic model, as it only requires the surface tilt angle as input, see \link{Koronakis}.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
V. Badescu, 3D isotropic approximation for solar diffuse irradiance on tilted surfaces, \emph{Renewable Energy}, Volume 26, Issue 2, 2002, Pages 221-233, ISSN 0960-1481, \url{http://dx.doi.org/10.1016/S0960-1481(01)00123-9}.
}

\seealso{
\code{\link{Liu}, \link{Koronakis}, \link{Tian}}
}
\examples{
data <- data.frame(s = 40)
Rd <- Badescu(data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Transposition}
\keyword{1st generation}
