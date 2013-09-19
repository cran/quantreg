\name{summary.crq}
\alias{summary.crqs}
\alias{summary.crq}
\alias{print.summary.crq}
\alias{print.summary.crqs}
\alias{plot.summary.crqs}
\title{
Summary methods for Censored Quantile Regression
}
\description{
Returns a summary object for a censored quantile regression fit.  A null value
will be returned if printing is invoked.
}
\usage{
\method{summary}{crq}(object, taus = 1:4/5, alpha = .05, se="boot", covariance=TRUE,  ...)
\method{print}{summary.crq}(x, digits = max(5, .Options$digits - 2), ...)
\method{print}{summary.crqs}(x,  ...)
\method{plot}{summary.crqs}(x, nrow = 3, ncol = 3, CoxPHit = NULL,  ...)
}
\arguments{
  \item{object}{
    An object of class \code{"crq"}  produced by a call to \code{crq()}.
  }
  \item{x}{
    An object of class \code{"crq"}  produced by a call to \code{crq()}.
  }
  \item{se}{
    specifies the method used to compute standard standard errors. but
    the only available method (so far) is "boot".  Further arguments to
    \code{\link{boot.crq}} and \code{\link{boot.rq}}  can be passed with 
    the \dots argument. 
  }
  \item{covariance}{
    logical flag to indicate whether the full covariance matrix of the 
    estimated parameters should be returned. 
  }
  \item{nrow}{Number of rows of the plot layout.}
  \item{ncol}{Number of columns of the plot layout.}
  \item{taus}{Quantiles to be summarized.}
  \item{alpha}{Confidence level for summary intervals.}
  \item{digits}{Number of digits to be printed in summary display.}
  \item{CoxPHit}{An object of class coxph produced by \code{coxph}.}
  \item{...}{
    Optional arguments to summary, e.g. to specify bootstrap methods
    sample sizes, etc.  see \code{\link{boot.rq}} and \code{\link{boot.crq}}
   }
}
\value{
  For method "FP" an object of class \code{summary.crqFP}  is returned 
  with the following components:
\item{coefficients}{
  a p by 4 matrix consisting of the coefficients, their estimated standard
  errors, their t-statistics, and their associated p-values.
}
\item{cov}{
  the estimated covariance matrix for the coefficients in the model,
  provided that \code{cov=TRUE} in the called sequence.
}
\item{rdf}{
  the residual degrees of freedom
}
\item{tau}{
  the quantile estimated
}
\item{cov}{
  the estimated covariance matrix for the coefficients in the model,
  provided that \code{cov=TRUE} in the called sequence.
}
  For method "PH" an object of class \code{summary.crq}  is returned 
  with the following components:
\item{coefficients}{
  a list of  p by 6 matrix consisting of the coefficients, upper and lower bounds
  for a (1-alpha) level confidence interval, their estimated standard
  errors, their t-statistics, and their associated p-values, one component for each
  element of the specified \code{taus} vector.
}
}
\details{
For the Fitzenberger-Powell method the resampling strategy used by the
\code{se = "boot"} method is based on the Bilias, Chen and Ying (2000) 
proposal.  For the Portnoy and Peng-Huang methods the bootstrapping 
is by default actually based on  a delete-d jackknife, as described in
Portnoy (2013), but resampling xy pairs using either conventional multinomial 
resampling or using exponential weighting as in Bose and Chatterjee (2003)
can be used by specifying the \code{bmethod} argument.  Note that the default
number of replications is set at \eqn{R = 100} a value that is obviously too small for
most applications.  This is done merely to speed up the examples in the 
documentation and facilitate testing. Larger, more appropriate values of \eqn{R}
can be passed to the bootstrapping functions via the \code{\ldots} argument
of the \code{summary} method.  It is important to recognize that when some
of the bootstrap replications are NA they are simply ignored in the computation
of the confidence bands and standard errors as currently reported.  The number
of these NAs is returned as part of the \code{summary.crq} object, and 
when printed is also reported. 
}

\references{
  Bose, A. and S. Chatterjee, (2003) Generalized bootstrap for estimators 
	of minimizers of convex functions, \emph{J. Stat. Planning and Inf}, 117,
	225-239.

  Bilias, Y. Chen, S. and Z. Ying, (2000) Simple resampling methods for censored
	quantile regression, \emph{J. of Econometrics}, 99, 373-386.

  Portnoy, S. (2013) The Jackknife's Edge:  Inference for Censored Quantile Regression,
        \emph{CSDA}, forthcoming.

}
\seealso{
  \code{\link{crq}}, \code{\link{QTECox}}
}
\keyword{regression}
