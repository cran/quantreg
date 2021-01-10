\name{boot.rq.pwxy}
\alias{boot.rq.pwxy}
\title{
    Preprocessing weighted bootstrap method
}
\description{
    Bootstrap method exploiting preprocessing strategy to reduce 
    computation time for large problem.  In contrast to 
    \code{\link{boot.rq.pxy}}  which uses the classical multinomial 
    sampling scheme and is coded in R, this uses the exponentially
    weighted bootstrap scheme and is coded in fortran and consequently
    is considerably faster in larger problems.
}
\usage{
boot.rq.pwxy(x, y, tau, coef, R = 200, m0 = NULL, eps = 1e-06, ...)
}
\arguments{
  \item{x}{
      Design matrix
}
  \item{y}{
      response vector
}
  \item{tau}{
      quantile of interest
}
  \item{coef}{
      point estimate of fitted object
}
  \item{R}{
      the number of bootstrap replications desired.
}
  \item{m0}{
      constant to determine initial sample size, defaults to sqrt(n*p)
      but could use some further tuning...
}
  \item{eps}{
      tolerance for convergence of fitting algorithm
}
  \item{...}{
      other parameters not yet envisaged.
}
}
\details{
    The fortran implementation is quite similar to the R code for 
    \code{\link{boot.rq.pxy}} except that there is no multinomial sampling.
    Instead \code{rexp(n)} weights are used.
}
\value{
  returns a list with elements: 
  \enumerate{
      \item{coefficients}{a matrix of dimension ncol(x) by R}
      \item{nit} {a 5 by m matrix of iteration counts}
      \item{info} {an m-vector of convergence flags}
  }
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
\code{\link{boot.rq.pxy}}
}
\keyword{bootstrap} 
