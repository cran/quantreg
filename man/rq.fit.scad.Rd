\name{rq.fit.scad}
\alias{rq.fit.scad}
\title{
SCADPenalized Quantile Regression 
}
\description{
  The fitting method implements the smoothly clipped absolute deviation
  penalty of Fan and Li for fitting quantile regression models.  
  When  the argument \code{lambda}
  is a scalar the penalty function is the scad modified l1
  norm of the last (p-1) coefficients, under the presumption that the
  first coefficient is an intercept parameter that should not be subject
  to the penalty.  When \code{lambda} is a vector it should have length
  equal the column dimension of the matrix \code{x} and then defines a
  coordinatewise specific vector of scad penalty parameters.  In this
  case \code{lambda} entries of zero indicate covariates that are not
  penalized.  There should be a sparse version of this, but isn't (yet).
}
\usage{
rq.fit.scad(x, y, tau = 0.5, alpha = 3.2, lambda = 1, start="rq", 
	beta = .9995, eps = 1e-06)
}
\arguments{
\item{x}{
  the design matrix
}
\item{y}{
  the response variable
}
\item{tau}{
  the quantile desired, defaults to 0.5.
}
\item{alpha}{
  tuning parameter of the scad penalty.
}
\item{lambda}{
  the value of the penalty parameter that determines how much shrinkage is done.
  This should be either a scalar, or a vector of length equal to the column dimension
  of the \code{x} matrix.

}
\item{start}{
  starting method, default method 'rq' uses the unconstrained rq estimate, while
  method 'lasso' uses the corresponding lasso estimate with the specified lambda.
}
\item{beta}{
  step length parameter for Frisch-Newton method.
}
\item{eps}{
 tolerance parameter for convergence. 
}
}
\value{
  Returns a list with a coefficient, residual, tau and lambda components.  
  When called from \code{"rq"} as intended the returned object
  has class "scadrqs".  
}
\details{The algorithm is an adaptation of the "difference convex algorithm"
described in Wu and Liu (2008).  It solves a sequence of (convex) QR problems 
to approximate solutions of the (non-convex) scad problem.} 
\references{
Wu, Y. and Y. Liu (2008) Variable  Selection in Quantile Regression, \emph{Statistica
Sinica}, to appear.
}
\author{R. Koenker}
\seealso{
\code{\link{rq}}}
\examples{
n <- 60
p <- 7
rho <- .5
beta <- c(3,1.5,0,2,0,0,0)
R <- matrix(0,p,p)
for(i in 1:p){
        for(j in 1:p){
                R[i,j] <- rho^abs(i-j)
                }
        }
set.seed(1234)
x <- matrix(rnorm(n*p),n,p) \%*\% t(chol(R))
y <- x \%*\% beta + rnorm(n)

f <- rq(y ~ x, method="scad",lambda = 30)
g <- rq(y ~ x, method="scad", start = "lasso", lambda = 30)
}
\keyword{regression}
