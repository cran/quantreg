\name{rrs.test}
\alias{rrs.test}
\title{
Quantile Regression Rankscore Test
}
\description{
  Function to compute regression rankscore test of a linear hypothesis
  based on the dual quantile regression process.  A test of  the 
  hypothesis, 
  is carried out by estimating the restricted model and constructing
  a test based on the dual process under the restricted model.  The
  details of the test are described in GJKP(1993).  The test has a
  Rao-score, Lagrange-multiplier interpretation since in effect it
  is based on the value of the gradient of unrestricted quantile regression
  problem evaluated under the null.  This function will eventually be
  superseded by a more general \code{anova()} method for \code{rq}.
}
\usage{
rrs.test(x0, x1, y, v, score="wilcoxon")
}
\arguments{
\item{x0}{
  the matrix of maintained regressors, a column of ones is
  appended automatically.
}
\item{x1}{
  matrix of covariates under test.
}
\item{y}{
  response variable, may be omitted if \code{v} is provided.
}
\item{v}{
  object of class \code{"rq.process"} generated e.g. by
  \code{rq(y ~ x0, tau=-1)}
}
\item{score}{
  Score function for test (see \code{\link{ranks}})
}
}
\value{
  Test statistic \code{sn} is asymptotically Chi-squared with rank(X1) dfs.
  The vector of ranks is also returned as component \code{rank}.
}
\details{
  See GJKP(1993)
}
\references{
  [1] Gutenbrunner, C., Jureckova, J., Koenker, R. and
  Portnoy, S. (1993)  Tests of linear hypotheses based on
  regression rank scores.   \emph{Journal of Nonparametric
  Statistics}, (2), 307-331.


  [2] Koenker, R. W. and d'Orey (1994).  Remark on Alg. AS
  229:  Computing dual regression quantiles and regression
  rank scores.  \emph{Applied Statistics}, \bold{43}, 410-414.
}
\seealso{
  \code{\link{rq}}, \code{\link{ranks}}
}
\examples{
# Test that covariates 2 and 3 belong in stackloss model using Wilcoxon scores.
data(stackloss)
rrs.test(stack.x[,1], stack.x[,2:3], stack.loss)
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
