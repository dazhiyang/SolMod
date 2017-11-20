\name{load.Perez.par}
\alias{load.Perez.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load Perez model parameters
}
\description{
This function loads Perez model parameters \eqn{F_1}, \eqn{F_2}, as well as ranges of the sky's clearness \eqn{\varepsilon} into workspace.
}
\usage{
load.Perez.par(model = "Perez90")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Specification of the version of Perez model. Choices are "\code{Perez90}", "\code{Perez87a}" or "\code{Perez87b}". Default to "\code{Perez90}".}
}
\details{
Three sets of Perez model coefficients are stored in this function. The model coefficients for "\code{Perez90}" version were reported by Perez et al. (1990). The other two sets of coefficients were reported by Perez et al. (1987). "\code{Perez87a}" corresponds to the set of coefficients for the 25 degrees circumsolar model; "\code{Perez87b}" corresponds to the set of coefficients for the point source circumsolar model.
}
\value{
The function returns a list containing the following.
  \item{F.mat}{8 by 6 matrix containing the Perez model parameters.}
  \item{epsilon.mat}{8 by 2 matrix containing the partitions of \eqn{\varepsilon}.}
}
\references{
Richard Perez, Pierre Ineichen, Robert Seals, Joseph Michalsky, Ronald Stewart, Modeling daylight availability and irradiance components from direct and global irradiance, \emph{Solar Energy}, Volume 44, Issue 5, 1990, Pages 271-289, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/0038-092X(90)90055-H}.\cr
Richard Perez, Robert Seals, Pierre Ineichen, Ronald Stewart, David Menicucci, A new simplified version of the Perez diffuse irradiance model for tilted surfaces, \emph{Solar Energy}, Volume 39, Issue 3, 1987, Pages 221-231, ISSN 0038-092X, \url{http://dx.doi.org/10.1016/S0038-092X(87)80031-2}.
}

\section{Warning}{The \eqn{\varepsilon} partitions in Perez et al. (1990) and Perez et al. (1987) are different. In addition, the \eqn{\varepsilon} calculation is also different: in Perez et al. (1990) it is given by: \deqn{\varepsilon = \frac{(D_h+I)/D_h + 1.041 z^3}{1+1.041 z^3},} while in Perez et al. (1987) it is given by: \deqn{\varepsilon = (D_h+I)/D_h.}}

\seealso{
\code{\link{Perez1}, \link{Perez2}, \link{Perez3}, \link{Perez4}, \link{NREL}}
}
\examples{
#load Perez90 parameters
tmp <- load.Perez.par(model = "Perez90")
tmp

#load Perez87b parameters
tmp <- load.Perez.par(model = "Perez87b")
tmp
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

