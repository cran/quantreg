\name{rqs.fit}
\alias{rqs.fit}
\title{Function to fit multiple response quantile regression models}
\usage{
rqs.fit(x, y, tau=0.5, tol = 0.0001)
}
\arguments{
\item{x}{
  the design matrix an n by p matrix.
}
\item{y}{
  the response variable as a n by m matrix
}
\item{tau}{
  the quantile desired, if tau lies outside (0,1) 
}
\item{tol}{
  tolerance parameter for Barrodale and Roberts exterior point method.
 }
}

\description{Function intended for multiple response quantile regression
called from \code{boot.rq} for wild bootstrap option.}
\seealso{  \code{\link{boot.rq}} }

\keyword{ regression }
