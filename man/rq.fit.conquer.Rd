\name{rq.fit.conquer}
\alias{rq.fit.conquer}
\title{Optional Fitting Method for Quantile Regression}
\description{
    This fitting method provides a link to the gradient descent 
    for convolution smoothed quantile regression problem implemented
    in the \pkg{conquer} package of He et al (2020).}
\usage{
rq.fit.conquer (x, y, tau=0.5, kernel = c("Gaussian", "uniform",
    "parabolic", "triangular"), h = 0, standardize = TRUE, tol = 1e-04,
    iteMax = 5000, ci = FALSE, alpha = 0.05, B = 1000)
}
\arguments{
  \item{x}{design matrix usually supplied via rq(), expected to
 have a intercept as the first column }
  \item{y}{ response vector usually supplied via rq() }
  \item{tau}{ quantile of interest }
  \item{kernel}{A character string specifying the choice of
          kernel function. Default is "Gaussian". Other choices are
          "uniform", "parabolic" or "triangular".} 
  \item{h}{The bandwidth parameter for kernel smoothing of the QR
      objective function.  Default is max{((log(n) + p) / n)^0.4, 0.05}. 
      The default is used if the input value is less than 0.05.} 
  \item{standardize}{A logical flag, when TRUE, the design matrix will be
          standardized so that each column has mean zero and standard
          deviation one. }
  \item{tol}{Tolerance level of the gradient descent
          algorithm. The gradient descent algorithm terminates when the
          maximal entry of the gradient is less than "tol". Default is
          1e-05.}
  \item{iteMax}{Maximum number of iterations. Default is 5000.} 
  \item{ci}{A logical flag. Default is FALSE. If "ci =
          TRUE", then three types of confidence intervals (percentile,
          pivotal and normal) will be constructed via multiplier
          bootstrap.  This option is subsumed in normal use by the
	  \code{summary.rq} functionality.}
  \item{alpha}{Nominal level for confidence intervals, may be passed
	via the call to \code{summary}}
  \item{B}{Number of bootstrap replications.  May be passed via summary.}
}
\details{
    See documentation in the \pkg{conquer} package.
}
\value{
Returns an object of class "rq".
}
\references{ 
Xuming He and Xiaoou Pan and Kean Ming Tan and Wen-Xin Zhou, (2020)
conquer: Convolution-Type Smoothed Quantile Regression,
\url{https://CRAN.R-project.org/package=conquer}}

\seealso{\code{\link{rq}}}

\keyword{regression}
