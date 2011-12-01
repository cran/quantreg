\name{rq.wfit}
\alias{rq.wfit}
\title{Function to choose method for Weighted Quantile Regression  }
\description{ Weight the data and then call the chosen fitting algorithm.  }
\usage{
rq.wfit(x, y, tau=0.5, weights, method="br", ...)
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
\item{weights}{ 
  weights used in the fitting
} 
\item{method}{ 
  method of computation:  "br" is Barrodale and Roberts exterior point 
  "fn" is the Frisch-Newton interior point method. 
} 
\item{...}{ 
  Optional arguments passed to fitting routine.
} 
}
 
\seealso{  \code{\link{rq}} \code{\link{rq.fit.br}} \code{\link{rq.fit.fnb}}}

\keyword{ regression }
