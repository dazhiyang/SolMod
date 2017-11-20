\name{Olmo2}
\alias{Olmo2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Olmo2 transposition model
}
\description{
Olmo's transposition model, modified by Evseev and Kudish (2009).
}
\usage{
Olmo2(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} containing zenith angle \eqn{z}, incidence angle \eqn{\theta}, GHI \eqn{G_h}, extraterrestrial GHI \eqn{I_{oh}} and albedo \eqn{\rho}.}
}
\details{
The original Olmo model, \link{Olmo1}, often performs poorly. Evseev and Kudish (2009) thus modified the model by considering sky conditions. Instead of modeling \eqn{G_c}, as in the original model, the modified version models \eqn{D_c}.
}
\value{
Returns the diffuse transposition factor \eqn{R_d}.
}
\references{
F.J Olmo, J Vida, I Foyo, Y Castro-Diez, L Alados-Arboledas, Prediction of global irradiance on inclined surfaces from horizontal global irradiance, \emph{Energy}, Volume 24, Issue 8, August 1999, Pages 689-704, ISSN 0360-5442, \url{http://dx.doi.org/10.1016/S0360-5442(99)00025-0}. \cr
Efim G. Evseev, Avraham I. Kudish, An assessment of a revised Olmo et al. model to predict solar global radiation on a tilted surface at Beer Sheva, Israel, \emph{Renewable Energy}, Volume 34, Issue 1, January 2009, Pages 112-119, ISSN 0960-1481, \url{http://dx.doi.org/10.1016/j.renene.2008.04.012}.
}


\seealso{
\code{\link{NREL}, \link{Olmo1}}
}
\examples{
data("NREL")
NREL$theta = NREL$theta1
NREL$s = 40
Rd <- Olmo2(NREL)

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

