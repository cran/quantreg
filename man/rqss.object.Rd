\name{rqss.object}
\alias{rqss.object}
\alias{logLik.rqss}
\alias{AIC.rqss}
\alias{fitted.rqss}
\alias{resid.rqss}
\alias{print.rqss}
\title{RQSS Objects and Summarization Thereof}
\description{
Functions to reveal the inner meaning of objects created by \code{rqss} fitting.
}
\usage{
\method{logLik}{rqss}(object, ...)
\method{AIC}{rqss}(object, ..., k=2)
\method{print}{rqss}(x, ...)
\method{resid}{rqss}(object, ...)
\method{fitted}{rqss}(object, ...)
}
\arguments{
  \item{object}{an object returned from \code{rqss} fitting, describing
	an additive model estimating a conditional quantile function. 
	See \code{\link{qss}} for details on how to specify these terms.}
  \item{x}{an rqss object, as above.}
  \item{k}{a constant factor governing the weight attached to the penalty
	term on effective degrees of freedom of the fit.  By default 
	k =2 corresponding to the Akaike version of the penalty, negative
	values indicate that the k should be set to log(n) as proposed
	by Schwarz (1978).} 
  \item{...}{additional arguments}
}
\details{ Total variation regularization for univariate and
    bivariate nonparametric quantile smoothing is described
    in Koenker, Ng and Portnoy (1994) and Koenker and Mizera(2003)
    respectively.  The additive model extension of this approach
    depends crucially on the sparse linear algebra implementation
    for R described in Koenker and Ng (2003).  Eventually, these
    functions should be expanded to provide an automated lambda
    selection procedure.}
\value{
    The function \code{summary.rqss} returns a list consisting of
    the following components:
    \item{fidelity}{Value of the quantile regression objective function.}
    \item{penalty}{A list consisting of the values of the total variation 
	smoothing penalty for each of additive components.}
    \item{edf}{Effective degrees of freedom of the fitted model, defined
	as the number of zero residuals of the fitted model,  Koenker
	Mizera (2003) for details.}
    \item{qssedfs}{A list of effective degrees of freedom for each of
	the additive components of the fitted model, defined as the
	number of non-zero elements of each penalty component of the
	residual vector.}
    \item{lamdas}{A list of the lambdas specified for each of the additive
	components of the model.}

}
\references{
  [1] Koenker, R. and S. Portnoy (1997)
  The Gaussian Hare and the Laplacean
  Tortoise:  Computability of Squared-error vs Absolute Error Estimators,
  (with discussion).
  \emph{Statistical Science} \bold{12}, 279--300.

  [2] Koenker, R., P. Ng and S. Portnoy, (1994)
  Quantile Smoothing Splines;
  \emph{Biometrika} \bold{81}, 673--680.

  [3] Koenker, R. and I. Mizera, (2003)
  Penalized Triograms: Total Variation Regularization for Bivariate Smoothing;
  \emph{JRSS(B)} \bold{66}, 145--163.

  [4] Koenker, R. and P. Ng (2003)
  SparseM:  A Sparse Linear Algebra Package for R,
  \emph{J. Stat. Software}.
}
\author{ Roger Koenker }
\seealso{ \code{\link{plot.rqss}}
}
\examples{
n <- 200
x <- sort(rchisq(n,4))
z <- x + rnorm(n)
y <- log(x)+ .1*(log(x))^2 + log(x)*rnorm(n)/4 + z
plot(x, y-z)
f.N  <- rqss(y ~ qss(x, constraint= "N") + z)
f.I  <- rqss(y ~ qss(x, constraint= "I") + z)
f.CI <- rqss(y ~ qss(x, constraint= "CI") + z)

lines(x[-1], f.N $coef[1] + f.N $coef[-(1:2)])
lines(x[-1], f.I $coef[1] + f.I $coef[-(1:2)], col="blue")
lines(x[-1], f.CI$coef[1] + f.CI$coef[-(1:2)], col="red")

## A bivariate example
data(CobarOre)
fCO <- rqss(z ~ qss(cbind(x,y), lambda= .08), data=CobarOre)
plot(fCO)
}
\keyword{regression}
\keyword{smooth}
\keyword{robust}
