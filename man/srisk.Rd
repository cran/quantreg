\name{srisk}
\alias{srisk}
\title{ Markowitz (Mean-Variance) Portfolio Optimization}
\description{
This function estimates optimal mean-variance portfolio weights from a matrix
of historical or simulated asset returns.  
}
\usage{
srisk(x, mu = 0.07, lambda = 1e+08, alpha = 0.1, eps = 1e-04)
}
\arguments{
  \item{x}{ Matrix of asset returns }
  \item{mu}{Required mean rate of return for the portfolio }
  \item{lambda}{Lagrange multiplier associated with mean return constraint}
  \item{alpha}{Choquet risk parameter, unimplemented }
  \item{eps}{ tolerance parameter for mean return constraint}
}
\details{
The portfolio weights are estimated by solving a constrained least squares problem.
}
\value{
  \item{pihat}{Optimal portfolio weights}
  \item{muhat }{Mean return in sample}
  \item{sighat }{Standard deviation of returns in sample}
}
\author{ R. Koenker }

\seealso{ \code{\link{qrisk}}}

\keyword{ regression }
