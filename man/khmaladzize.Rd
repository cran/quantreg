\name{khmaladzize}
\alias{khmaladzize}
\title{ Function to compute Khmaladze Transformation}
\description{
Function to compute the recursive least squares transformation
of the quantile regression process for the rq.test.khmal procedure.
}
\usage{
khmaladzize(tau, atau, Z, location.scale)
}
\arguments{
  \item{tau}{ quantiles specified in the fitted model}
  \item{atau}{ xbar'betahat(tau) at these quantiles}
  \item{Z}{ full rq process}
  \item{location.scale}{ if T do location-scale transformation,
			if F do location transformation}
}
\details{
Uses adaptive kernel density estimation akj() to estimate score functions.
}
\value{
Returns transformed Z process.
}
\references{ 
Koenker, Roger and Zhijie Xiao (2000), "Inference on the Quantile
     Regression Process'', unpublished. <URL:
     http://www.econ.uiuc.edu/~roger/research/inference/inference.html>
	}
\author{R. Koenker }
\seealso{  \code{\link{rq.test.khmal}}}

\keyword{htest}
