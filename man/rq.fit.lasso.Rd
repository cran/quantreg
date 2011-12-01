\name{rq.fit.lasso}
\alias{rq.fit.lasso}
\title{
Lasso Penalized Quantile Regression 
}
\description{
  The fitting method implements the lasso penalty of Tibshirani for
  fitting quantile regression models.  When  the argument \code{lambda}
  is a scalar the penalty function is the l1
  norm of the last (p-1) coefficients, under the presumption that the
  first coefficient is an intercept parameter that should not be subject
  to the penalty.  When \code{lambda} is a vector it should have length
  equal the column dimension of the matrix \code{x} and then defines a
  coordinatewise specific vector of lasso penalty parameters.  In this
  case \code{lambda} entries of zero indicate covariates that are not
  penalized.  There should be a sparse version of this, but isn't (yet).
}
\usage{
rq.fit.lasso(x, y, tau = 0.5, lambda = 1, beta = .9995, eps = 1e-06)
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
\item{lambda}{
  the value of the penalty parameter(s) that determine how much shrinkage is done.
  This should be either a scalar, or a vector of length equal to the column dimension
  of the \code{x} matrix.
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
  When called from \code{"rq"} (as intended) the returned object
  has class "lassorqs".  
}
\references{
Koenker, R. (2005 \emph{Quantile Regression}, CUP.
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

f <- rq(y ~ x, method="lasso",lambda = 30)
g <- rq(y ~ x, method="lasso",lambda = c(rep(0,4),rep(30,4)))
}
\keyword{regression}
