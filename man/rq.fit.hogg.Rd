\name{rq.fit.hogg}
\alias{rq.fit.hogg}
\title{weighted quantile regression fitting}
\description{
Function to estimate a regression mmodel by minimizing the weighted sum of several
quantile regression functions.  See Koenker(1984) for an asymptotic look at these
estimators.  This is a slightly generalized version of what Zou and Yuan (2008) call
composite quantile regression in that it permits weighting of the components of the
objective function and also allows linear constraints on the coefficients.
}
\usage{
rq.fit.hogg(x, y, taus = c(0.1, 0.3, 0.5), weights = c(0.7, 0.2, 0.1), 
	    R = NULL, r = NULL, beta = 0.99995, eps = 1e-06)
}
\arguments{
  \item{x}{design matrix}
  \item{y}{response vector }
  \item{taus}{quantiles getting positive weight}
  \item{weights}{weights assigned to the quantiles }
  \item{R}{optional matrix describing linear inequality constraints}
  \item{r}{optional vector describing linear inequality constraints}
  \item{beta}{step length parameter of the Frisch Newton Algorithm}
  \item{eps}{tolerance parameter for the Frisch Newton Algorithm}
}
\details{
Mimimizes a weighted sum of quantile regression objective functions using
the specified taus.  The model permits distinct intercept parameters at
each of the specified taus, but the slope parameters are constrained to
be the same for all taus.  This estimator was originally suggested to
the author by Bob Hogg in one of his famous blue notes of 1979.
The algorithm used to solve the resulting linear programming problems
is either the Frisch Newton algorithm described in Portnoy and Koenker (1997),
or the closely related algorithm described in Koenker and Ng(2002) that
handles linear inequality constraints.  See \code{\link{qrisk}} for illustration
of its use in portfolio allocation.
}
\value{
  \item{coefficients}{estimated coefficients of the model}
}
\references{
Zou, Hui and and Ming Yuan (2008)  Composite quantile regression and the
Oracle model selection theory, Annals of Statistics, 36, 1108--11120.

Koenker, R. (1984) A note on L-estimates for linear models, 
Stat. and Prob Letters, 2, 323-5.

Portnoy, S. and Koenker, R. (1997) The Gaussian Hare and the 
Laplacean Tortoise:  Computability of Squared-error vs Absolute Error Estimators, 
(with discussion).  Statistical Science, (1997) 12, 279-300.

Koenker, R. and Ng, P (2003) Inequality Constrained Quantile Regression, preprint.
 }
\author{ Roger Koenker }


\seealso{ \code{\link{qrisk}}}

\keyword{regression}
\keyword{ robust }
