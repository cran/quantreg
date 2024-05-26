\name{rq.fit.qfnb}
\alias{rq.fit.qfnb}
\title{
Quantile Regression Fitting via Interior Point Methods
}
\description{
This is a lower level routine called by \code{rq()} to compute quantile
regression parameters using the Frisch-Newton algorithm.  In contrast to
method "fn" it computes solutions for all the specified taus inside a
fortran loop.  See \code{\link{rq.fit.pfnb}} for further details on a more
efficient preprocessing method.
}
\usage{
rq.fit.qfnb(x, y, tau)
}
\arguments{
\item{x}{
The design matrix
}
\item{y}{
The response vector
}
\item{tau}{
The quantiles of interest, must lie in (0,1), be sorted and preferably equally
spaced.
}
}
\value{
  returns a list with elements consisting of 
      \item{coefficients}{a matrix of dimension ncol(x) by length(taus)
      }
      \item{nit}{a 3-vector of iteration counts}
      \item{info}{a convergence flag}
}
\details{
  The details of the Frisch-Newton algorithm are explained in Koenker and Portnoy (1997).
  The basic idea can be traced back to the log-barrier methods proposed by
  Frisch in the 1950's for linear programming.  But the current
  implementation is based on proposals by Mehrotra and others in the
  recent (explosive) literature on interior point methods for solving linear 
  programming problems.  This function replaces an earlier one \code{rq.fit.fn},
  which required the initial dual values to be feasible.  The current version allows the
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
