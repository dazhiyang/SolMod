\name{Koronakis}
\alias{Koronakis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Koronakis transposition model
}
\description{
Koronakis's transposition model (Koronakis, 1986)
}
\usage{
Koronakis(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}}
}
\details{
This model appears to be an isotropic transposition model as it only requires the surface tilt angle as input. However, the particular expression considers an anisotropic sky, e.g., for \eqn{s = 90^{\circ}}, the southern sky contributes \eqn{66.7\%} of the total sky-dome diffuse radiation. It is therefore a pseudo-isotropic model.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
Pericles S. Koronakis, On the choice of the angle of tilt for south facing solar collectors in the Athens basin area, \emph{Solar Energy}, Volume 36, Issue 3, 1986, Pages 217-225, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(86)90137-4}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Badescu}, \link{Liu}, \link{Tian}}
}
\examples{
data <- data.frame(s = 40)
Rd <- Koronakis(data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Transposition}
\keyword{1st generation}
