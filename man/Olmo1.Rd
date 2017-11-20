\name{Olmo1}
\alias{Olmo1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Olmo1 transposition model
}
\description{
Olmo's transposition model (Olmo et al., 1999).
}
\usage{
Olmo1(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing the tilt angle \eqn{s}, zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, DHI \eqn{D_h}, extraterrestrial GHI \eqn{I_{oh}} and albedo \eqn{\rho}.}
}
\details{
Olmo et al.'s model (Olmo et al., 1999) models the tilted global irradiance directly, instead of modeling the tilted diffuse irradiance.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
F.J Olmo, J Vida, I Foyo, Y Castro-Diez, L Alados-Arboledas, Prediction of global irradiance on inclined surfaces from horizontal global irradiance, \emph{Energy}, Volume 24, Issue 8, August 1999, Pages 689-704, ISSN 0360-5442, \url{http://dx.doi.org/10.1016/S0360-5442(99)00025-0}.
}


\seealso{
\code{\link{NREL}, \link{Olmo2}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Olmo1(NREL)

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


