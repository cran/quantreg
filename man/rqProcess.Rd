\name{rqProcess}
\alias{rqProcess}
\title{ Compute Standardized Quantile Regression Process }
\description{
Computes a standardize quantile regression process for the model
specified by the formula, on the partition of [0,1] specified by the
taus argument, and standardized according to the argument nullH. 
Intended for use in \code{\link{KhmaladzeTest}}.}
\usage{
rqProcess(formula, data, taus, nullH = "location", ...)
}
\arguments{
  \item{formula}{model formula  }
  \item{data}{data frame to be used to interpret formula }
  \item{taus}{ quantiles at which the process is to be evaluated, if any
	of the taus lie outside (0,1)  then the full process is computed 
	for all distinct solutions.}
  \item{nullH}{Null hypothesis to be used for standardization}
  \item{...}{optional arguments passed to \code{\link{summary.rq}}}
}
\details{
The process computes standardized estimates based on the 
hypothesis specified in the \code{nullH} argument.  
The Vhat component is rescaled by the Cholesky 
decomposition of the tau specific covariance matrix, the vhat component is
rescaled by the marginal standard errors.  The nature of the covariance
matrix used for the standardization is controlled arguments passed via
the \code{...} argument to \code{\link{summary.rq}}.  If the full
process is estimated then these covariance options aren't available
and only a simple iid-error form of the covariance matrix is used.
}
\value{
  \item{taus}{The points of evaluation of the process}
  \item{qtaus}{Values of xbar'betahat(taus)}
  \item{Vhat}{Joint parametric QR process}
  \item{vhat}{Marginal parametric QR processes}
}
\author{R. Koenker}
\seealso{\code{\link{KhmaladzeTest}}}
\keyword{regression}
