\name{rrs.test}
\alias{rrs.test}
\title{
Function to compute regression rankscore test of linear hypothesis
}
\description{
Tests the hypothesis that b_1 = 0 in the quantile regression model



	y = X_0 b_0 + X_1 b_1 + u 
}
\usage{
rrs.test(x0, x1, y, v, score="wilcoxon")
}
\arguments{
\item{x0}{
matrix of maintained regressors, a column of ones is appended automatically
}
\item{x1}{
matrix of regressors under test
}
\item{y}{
response variable may be omitted if v is provided
}
\item{v}{
regression quantile structure from rq(x0,y)
}
\item{score}{
score function for test (see ranks())
}}
\value{


\item{sn}{
Test statistic is asymptotically Chi-squared with rank(X1) degrees of freedom.
}
\item{rank}{
vector of ranks
}}
\section{Side Effects}{
}
\details{
}
\references{
[1] Gutenbrunner, C., J. Jureckova,  Koenker, R. and  Portnoy, S.(1993)
"Tests of Linear Hypotheses  based on Regression Rank Scores",
Journal of Nonparametric Statistics, (2), 307-331.



[2] Koenker, R.W. and d'Orey (1994).  "Remark on Alg. AS 229: Computing Dual
Regression Quantiles and Regression Rank Scores", Applied Statistics, 43, 410-414.



}
\seealso{
rq, ranks
}
\examples{
data(stackloss)
rrs.test(stack.x[,1:2],stack.x[,3],stack.loss)
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-2.
