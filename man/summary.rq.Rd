\name{summary.rq}
\alias{summary.rq}
\title{
Summary method for Quantile Regression
}
\description{
Returns a summary list for a quantile regression fit.  A null value
will be returned if printing is invoked.
}
\usage{
summary.rq(object, se="nid", covariance=TRUE, hs = TRUE,  ...)
}
\arguments{
  \item{object}{
    This is an object of class \code{"rq"} produced by a call to \code{rq()}.
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
    see \code{\link{boot.rq}}
   }
}
\value{
  a list is returned with the following components

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
  \code{se != "iid"}.
}
\item{J}{
  Outer product of gradient matrix returned if \code{cov=TRUE} and \code{se
    != "iid"}. The Huber sandwich is \code{cov = Hinv \%*\% J \%*\% Hinv}.
}}
\details{
When the default summary method is used, it tries to estimate a sandwich
form of the asymptotic covariance matrix and this involves estimating
the conditional density at each of the sample observations, negative
estimates can occur if there is crossing of the neighboring quantile
surfaces used to compute the difference quotient estimate.  If the
number of these is large relative to the sample size it is sometimes
an indication that some additional nonlinearity in the covariates
would be helpful, for instance quadratic effects.}
 

\references{
  Koenker, R. (2004) \emph{Quantile Regression}.
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
% Converted by Sd2Rd version 0.3-3.
