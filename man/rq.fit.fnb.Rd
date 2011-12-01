\name{rq.fit.fnb}
\alias{rq.fit.fnb}
\title{
Quantile Regression Fitting via Interior Point Methods
}
\description{
This is a lower level routine called by \code{rq()} to compute quantile
regression methods using the Frisch-Newton algorithm.
}
\usage{
rq.fit.fnb(x, y, tau=0.5, beta=0.99995, eps=1e-06)
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
  programming problems.  This function replaces an earlier one \code{rq.fit.fn},
  which required the initial dual values to be feasible.  This version allows the
  user to specify an infeasible starting point for the dual problem, that
  is one that may not satisfy the dual equality constraints.  It still 
  assumes that the starting value satisfies the upper and lower bounds.
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
