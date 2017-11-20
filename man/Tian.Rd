\name{Tian}
\alias{Tian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tian transposition model
}
\description{
Tian's transposition model (Tian et al., 2002)
}
\usage{
Tian(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}}
}
\details{
This is also a pseudo-isotropic model, see \link{Koronakis}.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Y.Q. Tian, R.J. Davies-Colley, P. Gong, B.W. Thorrold, Estimating solar radiation on slopes of arbitrary aspect, \emph{Agricultural and Forest Meteorology}, Volume 109, Issue 1, 2001, Pages 67-74, ISSN 0168-1923, \url{http://dx.doi.org/10.1016/S0168-1923(01)00245-3}.
}

\seealso{
\code{\link{Badescu}, \link{Koronakis}, \link{Liu}}
}
\examples{
data <- data.frame(s = 40)
Rd <- Tian(data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Transposition}
\keyword{1st generation}
