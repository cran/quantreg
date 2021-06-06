\name{LassoLambdaHat}
\alias{LassoLambdaHat}
\title{Lambda selection for QR lasso problems}
\description{
    Default procedure for selection of lambda in lasso constrained
    quantile regression as proposed by Belloni and Chernozhukov (2011)
}
\usage{
LassoLambdaHat(X, R = 1000, tau = 0.5, C = 1, alpha = 0.95)
}
\arguments{
  \item{X}{Design matrix}
  \item{R}{Number of replications}
  \item{tau}{quantile of interest}
  \item{C}{Cosmological constant}
  \item{alpha}{Interval threshold}
}
\value{
    vector of default lambda values of length p, the column dimension of X.
}
\details{
  As proposed by Belloni and Chernozhukov, a reasonable default lambda
  would be the upper quantile of the simulated values.  The procedure is based 
  on idea that a simulated gradient can be used as a pivotal statistic.
  Elements of the default vector are standardized by the respective standard deviations
  of the covariates. Note that the sqrt(tau(1-tau)) factor cancels in their (2.4) (2.6).
  In this formulation even the intercept is penalized.  If the lower limit of the
  simulated interval is desired one can specify \code{alpha = 0.05}.

} 
\references{
    Belloni, A. and  V. Chernozhukov. (2011) l1-penalized quantile regression 
    in high-dimensional sparse models. \emph{Annals of Statistics}, 39 82 - 130.
}
\examples{
n <- 200
p <- 10
x <- matrix(rnorm(n*p), n, p)
b <- c(1,1, rep(0, p-2))
y <- x \%*\% b + rnorm(n)
f <- rq(y ~ x, tau = 0.8, method = "lasso")
# See f$lambda to see the default lambda selection
}
