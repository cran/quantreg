\name{rq.fit.fnc}
\alias{rq.fit.fnc}
\title{
Quantile Regression Fitting via Interior Point Methods
}
\description{
This is a lower level routine called by \code{rq()} to compute quantile
regression methods using the Frisch-Newton algorithm.  It allows the
call to specify linear inequality constraints to which the fitted
coefficients will be subjected.  The constraints are assumed to be 
formulated as Rb >= r.
}
\usage{
rq.fit.fnc(x, y, R, r, tau=0.5, beta=0.9995, eps=1e-06)
}
\arguments{
\item{x}{
The design matrix
}
\item{y}{
The response vector
}
\item{R}{
The matrix describing the inequality constraints
}
\item{r}{
The right hand side vector of inequality constraints
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
  The details of the algorithm are explained in Koenker and Ng (2002).
  The basic idea can be traced back to the log-barrier methods proposed by
  Frisch in the 1950's for constrained optimization.  But the current
  implementation is based on proposals by Mehrotra and others in the
  recent (explosive) literature on interior point methods for solving linear 
  programming problems.  See \code{"rq"} helpfile for an example.
  It is an open research problem to provide an inference apparatus for
  inequality constrained quantile regression.
}
\references{
Koenker, R. and S. Portnoy (1997).
The Gaussian Hare and the Laplacian Tortoise:
Computability of squared-error vs. absolute-error estimators, with discussion,
\emph{Statistical Science}, \bold{12}, 279-300.

Koenker, R. and P. Ng(2005).
Inequality Constrained Quantile Regression, \emph{Sankya}, 418-440.
}
\seealso{
  \code{\link{rq}}, \code{\link{rq.fit.br}},
  \code{\link{rq.fit.pfn}}
}
\keyword{regression}
