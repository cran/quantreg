\name{rq.fit}
\alias{rq.fit}
\title{Function to choose method for Quantile Regression  }
\usage{
rq.fit(x, y, tau=0.5, method="br", ...)
}
\arguments{
\item{x}{
  the design matrix
}
\item{y}{
  the response variable
}
\item{tau}{
  the quantile desired, if tau lies outside (0,1) the whole process
  is estimated.
}
\item{method}{
  method of computation:  "br" is Barrodale and Roberts exterior point
  "fn" is the Frisch-Newton interior point method.
 }
\item{...}{
  Optional arguments passed to fitting routine.
 }
}

\description{Function to choose method for quantile regression}
\seealso{  \code{\link{rq}} \code{\link{rq.fit.br}} \code{\link{rq.fit.fnb}}}

\keyword{ regression }
