\name{lm.fit.recursive}
\alias{lm.fit.recursive}
\title{ Recursive Least Squares }
\description{
This function fits a linear model by recursive least squares.  It is
a utility routine for the rq.test.khmal function of the quantile regression
package.
}
\usage{
lm.fit.recursive(X, y, int=T)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{X}{ Design Matrix }
  \item{y}{ Response Variable}
  \item{int}{ if T then append intercept to X}
}
\value{
return p by n matrix of fitted parameters, where p. The
ith column gives the solution up to "time" i. 
}
\references{ A. Harvey, (1993) Time Series Models, MIT }
\author{ R. Koenker }
\seealso{khmaladzize }

\keyword{methods}
