\name{rq.fit.ppro}
\alias{rq.fit.ppro}
\title{
    Preprocessing fitting method for QR
}
\description{
    Preprocessing method for fitting quantile regression models that
    exploits the fact that adjacent tau's should have nearly the same
    sign vectors for residuals.
}
\usage{
rq.fit.ppro(x, y, tau, weights = NULL, Mm.factor = 0.8, eps = 1e-06, ...)
}
\arguments{
  \item{x}{
      Design matrix
}
  \item{y}{
      Response vector
}
  \item{tau}{
      quantile vector of interest 
}
  \item{weights}{
      case weights
}
  \item{Mm.factor}{
      constant determining initial sample size
}
  \item{eps}{
      Convergence tolerance
}
  \item{\dots}{
      Other arguments
}
}
\details{
    See references for further details.
}
\value{
    Returns a list with components:
    \item{coefficients}{Matrix of coefficient estimates}
    \item{residuals}{Matrix of residual estimates}
    \item{rho}{vector of objective function values}
    \item{weights}{vector of case weights}
}
\references{
     Chernozhukov, V.  I. Fernandez-Val and B. Melly,
	Fast Algorithms for the Quantile Regression Process, 2019,
	arXiv, 1909.05782,

     Portnoy, S.  and R. Koenker, The Gaussian Hare and the Laplacian
	Tortoise, Statistical Science, (1997) 279-300
}
\author{
    Blaise Melly and Roger Koenker
}
\seealso{
\code{\link{rq.fit.pfn}}, \code{\link{boot.rq.pxy}} 
}
\keyword{regression}
