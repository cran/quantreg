\name{trq.print}
\alias{trq.print}
\title{
Print a Trimmed Regression Quantile Summary
}
\description{
Prints out a summary of a trimmed quantile regression (using trq) including 
the standard error of the residuals, the overall F-statistic and t-statistics 
and standard errors for each coefficient.
}
\usage{
trq.print(trq.out, digits=4)
}
\arguments{
\item{trq.out}{
a list like the output from trq.
}
\item{digits}{
the number of significant digits that should be printed.
}}
\value{
an invisible list with the following components:

\item{summary}{
a vector containig the winsorized standard error of regression, 
the number of observations, the overall F-statistic for the regression, 
the degrees of freedom and the p-value for the F-statistic. 
When a1+a2=1, residual winsorized standard
error is the the square root of a Siddiqui-type estimate of the asymtotic
variance of the a1-quantile.
}
\item{coef.table}{
a matrix with columns containing the regression coefficients, 
their standard error, the t-statistic for testing if the coefficient 
is zero and the two sided p-value for the t statistic.
}}
\section{Side Effects}{
a summary of the regression (like the returned value) is printed out.
}
\seealso{
trq, lsfit, and ls.print for further reference.
}
\examples{
trq.print(trq(x,y))
}
\keyword{~keyword}
% Converted by Sd2Rd version 0.3-2.
