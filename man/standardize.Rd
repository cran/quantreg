\name{standardize}
\alias{standardize}
\title{ Function to standardize the quantile regression process}
\description{
Standardize the rq process prior to Khmaladze transformation.
}
\usage{
standardize(rqfit, location.scale=T)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{rqfit}{ Object produced by rq.test.khmal rq fitting}
  \item{location.scale}{ If T location-scale test, if F location test}
}
\details{
Compute standardized rq process.
}
\value{
Produces a list with several objects understood by rq.test.khmal
}
\references{ 
Koenker, Roger and Zhijie Xiao (2000), "Inference on the Quantile
     Regression Process'', unpublished. <URL:
     http://www.econ.uiuc.edu/~roger/research/inference/inference.html>
	}


\author{ R. Koenker }

\keyword{htest}
