\name{standardize}
\alias{standardize}
\title{ Function to standardize the quantile regression process}
\description{
Standardize the rq process prior to Khmaladze transformation.
}
\usage{
standardize(fit, nullH = "location-scale")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{fit}{ Object produced by rqProcess }
  \item{nullH}{ null hypothesis taking the value "location-scale" 
	or "location" }
}
\details{
Compute standardized rq process.
}
\value{
Produces a list with several objects understood by
\code{\link{khmaladze.test}}
}
\references{ 
Koenker, Roger and Zhijie Xiao (2000), "Inference on the Quantile
     Regression Process'', \textit{Econometrica},  81, 1583--1612. <URL:
     http://www.econ.uiuc.edu/~roger/research/inference/inference.html>
	}


\author{ R. Koenker }

\keyword{htest}
