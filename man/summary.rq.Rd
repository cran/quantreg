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
\method{summary}{rq}(object, se = NULL, covariance=FALSE, hs = TRUE,  U = NULL, gamma = 0.7, ...)
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
    are currently six available methods:  
    \enumerate{
      \item \code{"rank"} which produces confidence intervals for the
      estimated parameters by inverting a rank test as described in
      Koenker (1994).  This method involves solving a parametric linear
      programming problem, and for large sample sizes can be extremely
      slow, so by default it is only used when the sample size is less
      than 1000, see below.  The default option assumes that the errors are
      iid, while the option iid = FALSE implements a proposal of Koenker
      Machado (1999).  See the documentation for \code{rq.fit.br} for additional arguments.
      \item \code{"iid"} which presumes that the errors are iid and computes
      an estimate of the asymptotic covariance matrix as in KB(1978).
      
      \item \code{"nid"} which presumes local (in \code{tau})
      linearity (in \code{x}) of the
      the conditional quantile functions and computes a Huber
      sandwich estimate using a local estimate of the sparsity.
      
      \item \code{"ker"} which uses a kernel estimate of the sandwich
      as proposed by Powell(1990).

      \item \code{"boot"} which implements one of several possible bootstrapping
      alternatives for estimating standard errors including a variate of the wild
      bootstrap for clustered response.  See \code{\link{boot.rq}} for
      further details.  

      \item \code{"BLB"} which implements the bag of little bootstraps method
      proposed in Kleiner, et al (2014).  The sample size of the little bootstraps
      is controlled by the parameter \code{gamma}, see below.  At present only
      \code{bsmethod = "xy"} is sanction, and even that is experimental.  This
      option is intended for applications with very large n where other flavors
      of the bootstrap can be slow.
    }
    If \code{se = NULL} (the default)  and \code{covariance = FALSE}, and
    the sample size is less than 1001, then the "rank" method is used, 
    otherwise the "nid" method is used.
  }
  \item{covariance}{
    logical flag to indicate whether the full covariance matrix of the 
    estimated parameters should be returned. 
  }
  \item{hs}{
    Use Hall Sheather bandwidth for sparsity estimation
    If false revert to Bofinger bandwidth.
   }
   \item{U}{Resampling indices or gradient evaluations used for bootstrap,
       see \code{\link{boot.rq}}.}
   \item{gamma}{parameter controlling the effective sample size of the'bag
       of little bootstrap samples that will be \code{b = n^gamma} where
       \code{n} is the sample size of the original model.}
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
  provided that \code{cov=TRUE} in the called sequence.  This option
  is not available when se = "rank".
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
\item{B}{Matrix of bootstrap realizations.}
\item{U}{Matrix of bootstrap randomization draws.}
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
There are several options for alternative resampling methods.  When
\code{summary.rqs} is invoked, that is when \code{summary} is called
for a \code{rqs} object consisting of several \code{taus}, the \code{B}
components of the returned object can be used to construct a joint covariance
matrix for the full object.}
 

\references{
  Koenker, R. (2004) \emph{Quantile Regression}.

  Bilias, Y. Chen, S. and Z. Ying, Simple resampling methods for censored
	quantile regression, \emph{J. of Econometrics}, 99, 373-386.

  Kleiner, A., Talwalkar, A., Sarkar, P. and Jordan, M.I. (2014) A Scalable
	bootstrap for massive data, \emph{JRSS(B)}, 76, 795-816.

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
