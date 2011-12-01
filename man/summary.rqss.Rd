\name{summary.rqss}
\alias{summary.rqss}
\alias{print.summary.rqss}
\title{Summary of rqss fit}
\description{ Summary Method for a fitted rqss model.  }
\usage{\method{summary}{rqss}(object, cov = FALSE, ztol = 1e-5, ...) }
\arguments{
  \item{object}{an object returned from \code{rqss} fitting, describing
	an additive model estimating a conditional quantile function. 
	See \code{\link{qss}} for details on how to specify these terms.}
  \item{cov}{if TRUE return covariance matrix for the parametric components
	as \code{Vcov} and a list of covariance matrices for the nonparametric
	components as \code{Vqss}}
  \item{ztol}{Zero tolerance parameter used to determine the number of zero
	residuals indicating the estimated parametric dimension of the model,
	the so-called effective degrees of freedom.}
  \item{...}{additional arguments}
}
\details{ This is a highly experimental function intended to explore
inferential methods for rqss fitting.  The function is modeled after
\code{summary.gam} in Simon Wood's (2006)  \pkg{mgcv} package.  (Of course,
Simon should not be blamed for any deficiencies in the current implementation.
The basic idea is to condition on the lambda selection and construct
quasi-Bayesian credibility intervals based on normal approximation of
the "posterior," as computed using the Powell kernel estimate of the
usual quantile regression sandwich.  See \code{\link{summary.rq}} for
further details and references.
The function produces a conventional coefficient table with standard errors
t-statistics and p-values for the coefficients on the parametric part of the
model, and another table for additive nonparametric effects.  The latter
reports F statistics intended to evaluate the significance of these components
individually.  In addition the fidelity (value of the QR objective function
evaluated at the fitted model), the effective degrees of freedom, and the
sample size are reported.
}
\value{
    \item{coef}{Table of estimated coefficients and their standard errors,
	t-statistics, and p-values for the parametric components of the model}
    \item{qsstab}{Table of approximate F statistics, effective degrees of freedom
	and values of the penalty terms for each of the additive nonparametric 
	components of the model, and the lambda values assigned to each.} 
    \item{fidelity}{Value of the quantile regression objective function.}
    \item{tau}{Quantile of the estimated model}
    \item{formula}{formula of the estimated model}
    \item{edf}{Effective degrees of freedom of the fitted model, defined
	as the number of zero residuals of the fitted model,  see Koenker
	Mizera (2003) for details.}
    \item{n}{The sample size used to fit the model.}
    \item{Vcov}{Estimated covariance matrix of the fitted parametric component}
    \item{Vqss}{List of estimated covariance matrices of the fitted 
	nonparametric component}

}
\references{
  [1] Koenker, R., P. Ng and S. Portnoy, (1994)
  Quantile Smoothing Splines;
  \emph{Biometrika} \bold{81}, 673--680.

  [2] Koenker, R. and I. Mizera, (2003)
  Penalized Triograms: Total Variation Regularization for Bivariate Smoothing;
  \emph{JRSS(B)} \bold{66}, 145--163.

  [3] Wood, S. (2006) \emph{Generalized Additive Models}, Chapman-Hall.
}
\author{ Roger Koenker }
\seealso{ \code{\link{plot.rqss}} }
\examples{
n <- 200
x <- sort(rchisq(n,4))
z <- x + rnorm(n)
y <- log(x)+ .1*(log(x))^2 + log(x)*rnorm(n)/4 + z
f  <- rqss(y ~ qss(x) + z)
summary(f)
}
\keyword{regression}
\keyword{smooth}
\keyword{robust}
