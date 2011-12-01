\name{nlrq.control}
\alias{nlrq.control}
\title{ Set control parameters for nlrq }
\description{
Set algorithmic parameters for nlrq (nonlinear quantile regression function)
}
\usage{
nlrq.control(maxiter=100, k=2, InitialStepSize = 1, big=1e+20, eps=1e-07, beta=0.97)
}
\arguments{
  \item{maxiter}{maximum number of allowed iterations}
  \item{k}{the number of iterations of the Meketon algorithm to be calculated 
	in each step, usually 2 is reasonable, occasionally it may be helpful 
	to set k=1 }
  \item{InitialStepSize}{ Starting value in \code{optim} to determine the step
	length of iterations.  The default value of 1 is sometimes too optimistic.
	In such cases, the value 0 forces optim to just barely stick its toe in
	the water.}
  \item{big}{ a large scalar}
  \item{eps}{ tolerance for convergence of the algorithm }
  \item{beta}{ a shrinkage parameter which controls the recentering process 
	in the interior point algorithm. }
}
\seealso{  \code{\link{nlrq}} }

\keyword{ environment}
