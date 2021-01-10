\name{rq.fit.pfnb}
\alias{rq.fit.pfnb}
\title{
Quantile Regression Fitting via Interior Point Methods
}
\description{
This is a lower level routine called by \code{rq()} to compute quantile
regression parameters using the Frisch-Newton algorithm.  It uses a form
of preprocessing to accelerate the computations for situations in which
several taus are required for the same model specification.
}
\usage{
rq.fit.pfnb(x, y, tau, m0 = NULL, eps = 1e-06)
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
\item{m0}{
    An initial reduced sample size by default is set to be 
    \code{round((n * (log(p) + 1) )^(2/3)} this could be explored further
    to aid performance in extreme cases.
}
\item{eps}{A tolerance parameter intended to bound the confidence band entries
    away from zero.}
}
\value{
  returns a list with elements consisting of 
  \enumerate{
      \item{coefficients}{a matrix of dimension ncol(x) by length(taus)
      }
      \item{nit} {a 5 by m matrix of iteration counts: first two coordinates
	  of each column are the number of interior point iterations, the third is the 
	  number of observations in the final globbed sample size, and the last two 
	  are the number of fixups and bad-fixups respectively.  This is intended to
	  aid fine tuning of the initial sample size, m0.}
      \item{info} {an m-vector of convergence flags}
  }
}
\details{
  The details of the Frisch-Newton algorithm are explained in Koenker and Portnoy (1997),
  as is the preprocessing idea which is related to partial sorting and the algorithms
  such as \code{kuantile} for univariate quantiles that operate in time O(n).
  The preprocessing idea of exploiting nearby quantile solutions to accelerate
  estimation of adjacent quantiles is proposed in Chernozhukov et al (2020).
  This version calls a fortran version of the preprocessing algorithm that accepts
  multiple taus.  The preprocessing approach is also implemented for a single tau
  in \code{rq.fit.pfn} which may be regarded as a prototype for this function since
  it is written entirely in R and therefore is easier to experiment with.
}
\references{
Koenker, R. and S. Portnoy (1997).
The Gaussian Hare and the Laplacian Tortoise:
Computability of squared-error vs. absolute-error estimators, with discussion,
\emph{Statistical Science}, \bold{12}, 279-300.

Chernozhukov, V., I., Fernandez-Val, and Melly, B. (2020), `Fast algorithms for 
the quantile regression process', Empirical Economics, forthcoming.
}
\seealso{
  \code{\link{rq}}, \code{\link{rq.fit.br}},
  \code{\link{rq.fit.pfn}}
}
\keyword{regression}
