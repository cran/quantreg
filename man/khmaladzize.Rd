\name{khmaladzize}
\alias{khmaladzize}
\title{ Function to compute Khmaladze Transformation}
\description{
Function to compute the recursive least squares transformation
of the quantile regression process for the \code{\link{khmaladze.test}} procedure.
}
\usage{
khmaladzize(tau, atau, Z, nullH)
}
\arguments{
  \item{tau}{ quantiles specified in the fitted model}
  \item{atau}{ xbar'betahat(tau) at these quantiles}
  \item{Z}{ full rq process}
  \item{nullH}{ either "location-scale" 
			or "location" null hypothesis}
}
\details{
Uses adaptive kernel density estimation akj() to estimate score functions.
}
\value{
Returns transformed Z process.
}
\references{ 
Koenker, Roger and Zhijie Xiao (2001), "Inference on the Quantile
     Regression Process'', \textit{Econometrica},  81, 1583--1612. <URL:
     http://www.econ.uiuc.edu/~roger/research/inference/inference.html>
	}
\author{R. Koenker }
\seealso{  \code{\link{khmaladze.test}}}

\keyword{htest}
