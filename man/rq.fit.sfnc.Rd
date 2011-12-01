\name{rq.fit.sfnc}
\alias{rq.fit.sfnc}
\title{Sparse Constrained Regression Quantile Fitting}
\description{
  Fit constrained regression quantiles using a sparse implementation of
  the Frisch-Newton Interior-point algorithm.
}
\usage{
rq.fit.sfnc(x, y, R, r, tau = 0.5,
            rhs = (1-tau)*c(t(x) \%*\% rep(1,length(y))),control)
}
\arguments{
  \item{x}{structure of the design matrix X stored in csr format}
  \item{y}{response vector}
  \item{R}{constraint matrix stored in csr format}
  \item{r}{right-hand-side of the constraint}
  \item{tau}{desired quantile}
  \item{rhs}{the right-hand-side of the dual problem; regular users
    shouldn't need to specify this.}
  \item{control}{control paramters for fitting see \code{sfn.control}}
}
\details{
  This is a sparse implementation of the Frisch-Newton algorithm for
  constrained quantile regression described in Koenker and Portnoy (1996).
  The sparse matrix linear algebra is implemented through the functions
  available in the R package \pkg{SparseM}.
}
\value{
  \item{coef}{Regression quantile coefficients}
  \item{ierr}{Error code for the internal Fortran routine \code{srqfn}:
    \describe{
      \item{1:}{ insufficient work space in call to \code{extract}}
      \item{3:}{ insufficient storage in iwork when calling ordmmd}
      \item{4:}{ insufficient storage in iwork when calling sfinit}
      \item{5:}{ nnzl > nnzlmax when calling sfinit}
      \item{6:}{ nsub > nsubmax when calling sfinit}
      \item{7:}{ insufficient work space in iwork when calling symfct}
      \item{8:}{ inconsistancy in input when calling symfct}
      \item{9:}{ tmpsiz > tmpmax when calling symfct; increase tmpmax}
      \item{10:}{ nonpositive diagonal encountered when calling blkfct}
      \item{11:}{ insufficient work storage in tmpvec when calling blkfct}
      \item{12:}{ insufficient work storage in iwork when calling blkfct}
      \item{13:}{ nnzd > nnzdmax in e,je when calling amub}
      \item{14:}{ nnzd > nnzdmax in g,jg when calling amub}
      \item{15:}{ nnzd > nnzdmax in h,jh when calling aplb}
      \item{15:}{ tiny diagonals replaced with Inf when calling blkfct}
    }
  }
  \item{it}{Iteration count}
  \item{time}{Amount of time used in the computation}
}
\references{
  Koenker, R and Ng, P. (2002).
  SparseM:  A Sparse Matrix Package for \R; \cr
  \url{http://www.econ.uiuc.edu/~roger/research}

  Koenker, R. and P. Ng(2005).
  Inequality Constrained Quantile Regression, \emph{Sankya}, 418-440.

}
\author{Pin Ng}
\seealso{
  \code{\link{rq.fit.sfn}} for the unconstrained version,
  \pkg{SparseM} for the underlying sparse matrix \R package.
}
\examples{
## An artificial example :
n <- 200
p <- 50
set.seed(17)
X <- rnorm(n*p)
X[abs(X) < 2.0] <- 0
X <- cbind(1,matrix(X,n,p))
y <- 0.5 * apply(X,1,sum) + rnorm(n) ## true beta = (0.5, 0.5, ...)
R <- rbind(diag(p+1), -diag(p+1))
r <- c(rep( 0, p+1), rep(-1, p+1))

sX <- as.matrix.csr(X)
sR <- as.matrix.csr(R)
try(rq.o <- rq.fit.sfnc(sX, y, sR, r)) #-> not enough tmp memory

(tmpmax <- floor(1e5 + exp(-12.1)*(sX@ia[p+1]-1)^2.35))
## now ok:
rq.o <- rq.fit.sfnc(sX, y, sR, r, control = list(tmpmax = tmpmax))
}
\keyword{regression}
