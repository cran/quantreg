\name{rq.test.khmal}
\alias{rq.test.khmal}
\title{
Estimates Quantile Regression Model and Test Statistics
}
\description{Estimates a model and produces the output necessary to test
the location and location-scale shift hypotheses. Returns an object of
the class "khmal".}

\usage{
rq.test.khmal( formula, data, taus=seq(0.2,0.8,by=0.002), 
location.scale = TRUE,  trim = c(0.25, 0.75) ) 
}
\arguments{
\item{formula}{
a symbolic description of the model to be fit. The details of model
specification are given below.
}
\item{data}{
dataframe containing the regressand and regressors. It should \bold{not} 
contain a column of ones.
}
\item{taus}{
vector of quantiles to be estimated. Default range is 0.2 to 0.8,
increasing by 0.002 units.
}
\item{location.scale}{
a logic value indicating whether the location-scale shift hypothesis
(default) or the location shift hypothesis (F) should be tested.
}
\item{trim}{
a vector indicating the lower and upper bound of the quantiles to
included in the computation of the test statistics (only, not
estimates). This might be required due to tail behavior.
}
}
\value{
a R structure is returned

\item{formula}{
a symbolic description of the fitted model.
}
\item{taus}{
vector with estimated quantiles.
}
\item{Jn}{
\eqn{X'X} of least squares.
}
\item{fit}{
Estimated coefficients for all quantiles.
}
\item{Hfit}{
Array with inverse quantile regression covariance matrices.
}
\item{vars}{
Variables names, useful for plotting functions.
}
\item{var.list}{
Numerical list of variables to be plotted. By default all variables
are plotted. A restricted set of variables can be specified by
providing a numerical vector indicating the desired variables. The
convention is that 1 corresponds to the intercept, 2 to the first 
independent variable entered in "formula" and so on. 
See example for further details.
}
\item{location.scale}{
a logic value indicating whether the location-scale shift hypothesis
(default) or the location shift hypothesis (F) was performed.
}
\item{b}{
Matrix with the coefficients from the least-squares regression of the
estimated quantile regression slopes on the intercept.
}
\item{J}{
Output of standardize().
}
\item{Vtilde}{
Output of khmaladzize() for the joint hypothesis testing.
}
\item{vtilde}{
Output of khmaladzize() for the individual hypotheses testing.
}
\item{Tvtilde}{
vtilde appropriately transformed to compute the individual test
statistics of the location-scale shift hypothesis.
}
\item{trim}{
a vector indicating the lower and upper bound of the quantiles
included in the computation of the test statistics (only, not
estimates).
}
\item{Kn}{
Joint location-scale shift test statistic. Only returned if
location.scale is TRUE.
}
\item{KHn}{
Individual location-scale shift test statistics. Only returned if
location.scale is TRUE.
}
\item{Tn}{
Joint location shift test statistic. Only returned if location.scale
is FALSE.
}
\item{THn}{
Individual location shift test statistics. Only returned if
location.scale is FALSE.
}
}

\examples{
data(barro)
fit.Lonly <- rq.test.khmal( y.net ~ lgdp2 + fse2 + gedy2 + Iy2 + gcony2, 
data = barro, location.scale = FALSE )
}
\keyword{htest}
\references{
Koenker, Roger and Zhijie Xiao (2000), "Inference on the Quantile
Regression Process'', unpublished.
\url{http://www.econ.uiuc.edu/~roger/research/inference/inference.html}
}
