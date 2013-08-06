\name{KhmaladzeTest}
\alias{KhmaladzeTest}
\alias{khmaladzize}
\alias{qdensity}
\title{ Tests of Location and Location Scale Shift Hypotheses for Linear Models}
\description{Tests of the hypothesis that a linear model specification
is of the location shift or location-scale shift form.  The tests are based 
on the Doob-Meyer Martingale transformation approach proposed by Khmaladze(1981)
for general goodness of fit problems, and adapted to quantile regression by 
Koenker and Xiao (2002).}

\usage{
KhmaladzeTest(formula, data = NULL, taus = -1, nullH = "location" ,  
	trim = c(0.05, 0.95), ...) 
}
\arguments{
\item{formula}{a formula specifying the model to fit by \code{\link{rqProcess}}}
\item{data}{a data frame within which to interpret the formula}
\item{taus}{An equally spaced grid of points on which to evaluate the 
quantile regression process, if any taus fall outside (0,1) then the full
process is computed.}
\item{nullH}{a character vector indicating whether the "location" shift hypothesis
(default) or the "location-scale" shift hypothesis  should be tested.  }
\item{trim}{ a vector indicating the lower and upper bound of the quantiles to
included in the computation of the test statistics (only, not
estimates).  }
\item{...}{other arguments to be passed to \code{\link{summary.rq}.}}
}
\value{
an object of class KhmaladzeTest  is returned containing:

\item{nullH}{ The form of the null hypothesis.}

\item{Tn}{
Joint test statistic of the hypothesis that all the slope
parameters of the model satisfy the hypothesis.
}
\item{THn}{
Vector of test statistics testing whether individual slope
parameters satisfy the null hypothesis.
}
}

\examples{
data(barro)
KhmaladzeTest( y.net ~ lgdp2 + fse2 + gedy2 + Iy2 + gcony2, 
		data = barro, taus = seq(.05,.95,by = .01))
}
\keyword{htest}
\references{
Khmaladze, E. (1981) ``Martingale Approach in the Theory of
Goodness-of-fit Tests,'' \emph{Theory of Prob. and its Apps}, 26,
240--257.

Koenker, Roger and Zhijie Xiao (2002), ``Inference on the Quantile
Regression Process'',  \emph{Econometrica},  81, 1583--1612.
\url{http://www.econ.uiuc.edu/~roger/research/inference/inference.html}
}
