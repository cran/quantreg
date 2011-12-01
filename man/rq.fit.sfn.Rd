\name{rq.fit.sfn}
\alias{rq.fit.sfn}
\alias{sfnMessage}
\title{Sparse Regression Quantile Fitting}
\description{
  Fit a quantile regression model using a sparse implementation of the
  Frisch-Newton interior-point algorithm.
}
\usage{
rq.fit.sfn(a, y, tau = 0.5, rhs = (1-tau)*c(t(a) \%*\% rep(1,length(y))), control)
}
\arguments{
  \item{a}{structure of the design matrix X stored in csr format}
  \item{y}{response vector}
  \item{tau}{desired quantile}
  \item{rhs}{the right-hand-side of the dual problem; regular users
    shouldn't need to specify this, but in special cases can be quite
    usefully altered to meet special needs.  See e.g. Section 6.8 of
    Koenker (2005).}
  \item{control}{control parameters for fitting routines:  see \code{sfn.control}}
}
\details{
  This is a sparse implementation of the Frisch-Newton algorithm for quantile
  regression described in Portnoy and Koenker (1997). The sparse matrix
  linear algebra is implemented through the functions available in the R
  package \pkg{SparseM}.
}
\value{
  \item{coef}{Regression quantile coefficients}
  \item{ierr}{Error code for the internal Fortran routine \code{srqfnc}:
    \describe{
      \item{1:}{ insufficient work space in call to \code{extract}}
      \item{2:}{ nnzd > nnzdmax}
      \item{3:}{ insufficient storage in iwork when calling ordmmd}
      \item{4:}{ insufficient storage in iwork when calling sfinit}
      \item{5:}{ nnzl > nnzlmax when calling sfinit}
      \item{6:}{ nsub > nsubmax when calling sfinit}
      \item{7:}{ insufficient work space in iwork when calling symfct}
      \item{8:}{ inconsistancy in input when calling symfct}
      \item{9:}{ tmpsiz > tmpmax when calling bfinit; increase tmpmax}
      \item{10:}{ nonpositive diagonal encountered, not positive definite}
      \item{11:}{ insufficient work storage in tmpvec when calling blkfct}
      \item{12:}{ insufficient work storage in iwork when calling blkfct}
      \item{17:}{ tiny diagonals replaced with Inf when calling blkfct}
    }
  }
  \item{it}{Iteration count}
  \item{time}{Amount of time used in the computation}
}
\references{
  Portnoy, S. and R. Koenker (1997) The Gaussian Hare and the Laplacean Tortoise: 
  Computability of Squared-error vs Absolute Error Estimators, (with discussion). 
  \emph{Statistical Science}, 12, 279-300.

  Koenker, R and Ng, P. (2003).  SparseM:  A Sparse Matrix Package for \R, 
  \emph{J. of Stat. Software}, 8, 1--9.

  Koenker, R. (2005) \emph{Quantile Regression}, Cambridge U. Press.
}
\author{Pin Ng}
\seealso{\code{rq.fit.sfnc} for the constrained version,
  \code{SparseM} for a sparse matrix package for \R
}
\examples{
## An artificial example :
n <- 200
p <- 50
set.seed(101)
X <- rnorm(n*p)
X[abs(X) < 2.0] <- 0
X <- cbind(1, matrix(X, n, p))
y <- 0.5 * apply(X,1,sum) + rnorm(n) ## true beta = (0.5, 0.5, ...)

sX <- as.matrix.csr(X)
try(rq.o <- rq.fit.sfn(sX, y)) #-> not enough tmp memory
(tmpmax <- floor(1e5 + exp(-12.1)*(sX@ia[p+1]-1)^2.35))
## now ok:
rq.o <- rq.fit.sfn(sX, y, control = list(tmpmax= tmpmax))
}
\keyword{regression}
