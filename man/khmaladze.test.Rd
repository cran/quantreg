\name{khmaladze.test}
\alias{khmaladze.test}
\title{ Tests of Location and Location Scale Hypothesis for Linear Models}
\description{Tests of the hypothesis that a linear model specification
is of the location and location-scale shift form.  The tests are based 
on the Doob-Meyer transformation approach proposed by Khmaladze(1981)
for general goodness of fit problems,
and adapted to quantile regression by Koenker and Xiao (2001).}

\usage{
khmaladze.test( fit, nullH = "location-scale" ,  trim = c(0.25, 0.75) ) 
}
\arguments{
\item{fit}{
an object produced by \code{\link{rqProcess}} containing
components describing the quantile regression process for
the model.
}
\item{nullH}{
a character vector indicating whether the "location-scale" shift hypothesis
(default) or the "location" shift hypothesis  should be tested.
}
\item{trim}{
a vector indicating the lower and upper bound of the quantiles to
included in the computation of the test statistics (only, not
estimates). This might be required due to tail behavior.
}
}
\value{
an object of class khmaladze  is returned containing:

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
fit <- rqProcess( y.net ~ lgdp2 + fse2 + gedy2 + Iy2 + gcony2, 
		data = barro, taus = seq(.1,.9,by = .05))
khmaladze.test(fit, nullH = "location")
}
\keyword{htest}
\references{
Khmaladze, E. (1981) ``Martingale Approach in the Theory of
Goodness-of-fit Tests,'' \textit{Theory of Prob. and its Apps}, 26,
240--257.

Koenker, Roger and Zhijie Xiao (2000), "Inference on the Quantile
Regression Process'',  \textit{Econometrica},  81, 1583--1612.
\url{http://www.econ.uiuc.edu/~roger/research/inference/inference.html}
}
