\name{rq.fit.fn}
\alias{rq.fit.fn}
\title{
Quantile Regression Fitting via Interior Point Methods
}
\description{
This is a lower level routine called by \code{rq()} to compute quantile
regression methods using the Frisch-Newton algorithm.
}
\usage{
rq.fit.fn(x, y, tau=0.5, int=F, beta=0.99995, eps=1e-06)
}
\arguments{
\item{x}{
The design matrix
}
\item{y}{
The response vector
}
\item{tau}{
The quantile of interest, must lie in (0,1)
}
\item{int}{
logical flag, if T an intercept should be appended to \code{x}, if not, not.
}
\item{beta}{
technical step length parameter -- alter at your own risk!
}
\item{eps}{
  tolerance parameter for convergence.  In cases of multiple optimal solutions
  there may be some discrepancy between solutions produced by method
  \code{"fn"} and method \code{"br"}.  This is due to the fact that
  \code{"fn"} tends to converge to a point near the centroid of the
  solution set, while \code{"br"} stops at a vertex of the set.  
}
}
\value{
  returns an object of class \code{"rq"}, which can be passed to
  \code{\link{summary.rq}} to obtain standard errors, etc.
}
\details{
  The details of the algorithm are explained in Koenker and Portnoy (1997).
  The basic idea can be traced back to the log-barrier methods proposed by
  Frisch in the 1950's for constrained optimization.  But the current
  implementation is based on proposals by Mehrotra and others in the
  recent (explosive) literature on interior point methods for solving linear 
  programming problems.  This version of the algorithm is designed for
  fairly large problems, for very large problems see \code{rq.fit.pfn}.
}
\references{
Koenker, R. and S. Portnoy (1997).
The Gaussian Hare and the Laplacian Tortoise:
Computability of squared-error vs. absolute-error estimators, with discussion,
\emph{Statistical Science}, \bold{12}, 279-300.
}
\seealso{
  \code{\link{rq}}, \code{\link{rq.fit.br}},
  \code{\link{rq.fit.pfn}}
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
