\name{akj}
\alias{akj}
\title{Density Estimation using Adaptive Kernel method}
\description{
  Univariate \emph{adaptive} kernel density estimation a la Silverman.
  As used by Portnoy and Koenker (1989).
}
\usage{
akj(x, z =, p =, h = -1, alpha = 0.5, kappa = 0.9, iker1 = 0)
}
\arguments{
  \item{x}{points used for centers of kernel assumed to be sorted.}
  \item{z}{points at which density is calculated; defaults to an
    equispaced sequence covering the range of x.}
  \item{p}{vector of probabilities associated with \code{x}s; defaults
    to 1/n for each x.}
  \item{h}{initial window size (overall); defaults to Silverman's normal
    reference.}
  \item{alpha}{a sensitivity parameter that determines the sensitivity of
    the local bandwidth to variations in the pilot density; defaults to .5.}
  \item{kappa}{constant multiplier for initial (default) window width}
  \item{iker1}{integer kernel indicator: 0 for normal kernel (default)
    while 1 for Cauchy kernel (\code{\link{dcauchy}}).}
}
\value{
  a \code{\link{list}} structure is with components
  \item{dens}{the vector of estimated density values \eqn{f(z)}}
  \item{psi}{a vector of \eqn{\psi=-f'/f} function values.}
  \item{score}{a vector of score \eqn{\psi' = (f'/f)^2 - f''/f} function
    values.}
  \item{h}{same as the input argument h}
}
\note{
  if the \code{score} function values are of interest, the Cauchy kernel
  may be preferable.
}
\references{
  Portnoy, S and R Koenker, (1989)
  Adaptive L Estimation of Linear Models;
  \emph{Annals of Statistics} \bold{17}, 362--81.

  Silverman, B. (1986)
  \emph{Density Estimation}, pp 100--104.
}
\examples{
 set.seed(1)
 x <- c(rnorm(600), 2 + 2*rnorm(400))
 xx <- seq(-5, 8, length=200)
 z <- akj(x, xx)
 plot(xx, z$dens, ylim=range(0,z$dens), type ="l", col=2)
 abline(h=0, col="gray", lty=3)
 plot(xx, z$psi, type ="l", col=2, main = expression(hat(psi(x))))
 plot(xx, z$score, type ="l", col=2,
      main = expression("score " * hat(psi) * "'" * (x)))

 if(require("nor1mix")) {
  m3 <- norMix(mu= c(-4, 0, 3), sig2 = c(1/3^2, 1, 2^2),
               w = c(.1,.5,.4))
  plot(m3, p.norm = FALSE)
  set.seed(11)
  x <- rnorMix(1000, m3)
  z2 <- akj(x, xx)
  lines(xx, z2$dens, col=2)
  z3 <- akj(x, xx, kappa = 0.5, alpha = 0.88)
  lines(xx, z3$dens, col=3)
 }
}
\keyword{smooth}
