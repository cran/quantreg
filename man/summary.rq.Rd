\name{summary.rq}
\alias{summary.rq}
\alias{summary.rqs}
\alias{summary.rcrqs}
\title{
Summary methods for Quantile Regression
}
\description{
Returns a summary list for a quantile regression fit.  A null value
will be returned if printing is invoked.
}
\usage{
\method{summary}{rq}(object, se = NULL, covariance=FALSE, hs = TRUE,  ...)
\method{summary}{rqs}(object, ...)
}
\arguments{
  \item{object}{
    This is an object of class \code{"rq"} or \code{"rqs"} produced by 
    a call to \code{rq()}, depending on whether one or more taus are
    specified.
  }
  \item{se}{
    specifies the method used to compute standard standard errors.  There
    are currently five available methods:  
    \enumerate{
      \item \code{"rank"} which produces confidence intervals for the
      estimated parameters by inverting a rank test as described in
      Koenker (1994).  The default option assumes that the errors are
      iid, while the option iid = FALSE implements the proposal of Koenker
      Machado (1999).  
      See the documentation for \code{rq.fit.br} for additional arguments.

      \item \code{"iid"} which presumes that the errors are iid and computes
      an estimate of the asymptotic covariance matrix as in KB(1978).
      
      \item \code{"nid"} which presumes local (in \code{tau})
      linearity (in \code{x}) of the
      the conditional quantile functions and computes a Huber
      sandwich estimate using a local estimate of the sparsity.
      
      \item \code{"ker"} which uses a kernel estimate of the sandwich
      as proposed by Powell(1990).

      \item \code{"boot"} which implements one of several possible bootstrapping
      alternatives for estimating standard errors.
    }
    If \code{se = NULL} (the default)  and \code{covariance = FALSE}
    then the "rank" method is used, otherwise the "nid" method is used.
  }
  \item{covariance}{
    logical flag to indicate whether the full covariance matrix of the 
    estimated parameters should be returned. 
  }
  \item{hs}{
    Use Hall Sheather bandwidth for sparsity estimation
    If false revert to Bofinger bandwidth.
   }
  \item{...}{
    Optional arguments to summary, e.g. bsmethod to use bootstrapping.
    see \code{\link{boot.rq}}.  When using the "rank" method for confidence
    intervals, which is the default method for sample sizes less than 1000,
    the type I error probability of the intervals can be controlled with the
    alpha parameter passed via "...",  thereby controlling the width of the
    intervals plotted by \code{plot.summary.rqs}. 
   }
}
\value{
  a list is returned with the following components, when \code{object}
  is of class \code{"rqs"} then there is a list of such lists.

\item{coefficients}{
  a p by 4 matrix consisting of the coefficients, their estimated standard
  errors, their t-statistics, and their associated p-values.
}
\item{cov}{
  the estimated covariance matrix for the coefficients in the model,
  provided that \code{cov=TRUE} in the called sequence.
}
\item{Hinv}{
  inverse of the estimated Hessian matrix returned if \code{cov=TRUE} and
  \code{se \%in\% c("nid","ker") }, note that for \code{se = "boot"} there
  is no way to split the estimated covariance matrix into its sandwich
  constituent parts.
}
\item{J}{
  Unscaled Outer product of gradient matrix returned if \code{cov=TRUE} and \code{se
    != "iid"}. The Huber sandwich is \code{cov = tau (1-tau) Hinv \%*\% J \%*\% Hinv}.
    as for the \code{Hinv} component, there is no \code{J} component when
    \code{se == "boot"}.  (Note that to make the Huber sandwich you need to add the 
    tau (1-tau) mayonnaise yourself.)
    }
}
\details{
When the default summary method is used, it tries to estimate a sandwich
form of the asymptotic covariance matrix and this involves estimating
the conditional density at each of the sample observations, negative
estimates can occur if there is crossing of the neighboring quantile
surfaces used to compute the difference quotient estimate.  
A warning message is issued when such negative estimates exist indicating
the number of occurrences -- if this number constitutes a large proportion
of the sample size, then it would be prudent to consider an alternative 
inference method like the bootstrap.
If the number of these is large relative to the sample size it is sometimes
an indication that some additional nonlinearity in the covariates
would be helpful, for instance quadratic effects.
Note that the default \code{se} method is rank, unless the sample size exceeds
1001, in which case the \code{rank} method is used.
}
 

\references{
  Koenker, R. (2004) \emph{Quantile Regression}.

  Bilias, Y. Chen, S. and Z. Ying, Simple resampling methods for censored
	quantile regression, \emph{J. of Econometrics}, 99, 373-386.
}
\seealso{
  \code{\link{rq}}
  \code{\link{bandwidth.rq}}
}
\examples{
data(stackloss)
y <- stack.loss
x <- stack.x
summary(rq(y ~ x, method="fn")) # Compute se's for fit using "nid" method.
summary(rq(y ~ x, ci=FALSE),se="ker")
# default "br" alg, and compute kernel method se's
}
\keyword{regression}
