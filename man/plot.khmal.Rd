\name{plot.khmal}
\alias{plot.khmal}
\title{
Plots Standardized and Khmaladzized Residual Processes 
}
\description{
The function makes 6 arrays of p plots based on the object of class
"khmal" created by `rq.test.khmal' of quantile regression results. The 6
arrays are: (i) estimated coefficients; (ii) regression of slopes on
the intercept; (iii) standardized residuals for the joint; (iv)
standardized residuals for the coef by coef; (v) khmaladzized
residuals for the joint, and (vi) khmaladzized residuals for the coef
by coef hypothesis testing.
} 

\usage{
plot.khmal(x, nrow= ceiling(length(x$var.list)/2), ncol= 2, plotn = 1:6, 
	bcolor="gray", \dots )
}
\arguments{
\item{x}{
output of `rq.test.khmal'. plot.khmal() \bold{requires} the output of 
`rq.test.khmal'.
}
\item{var.list}{
numerical list of variables to be plotted. By default all variables are 
plotted. A restricted set of variables can be specified by providing a 
numerical vector indicating the desired variables. The convention is 
that 1 corresponds to the intercept, 2 to the first independent variable 
entered in "formula" and so on. See example for further details.
}
\item{nrow}{
number of rows per page of plots. Automatically set by assuming that
the number of columns is 2.
}
\item{ncol}{
number of plots per page of plots. Default 2.
}
\item{plotn}{
a numerical vector indicating which array of plots will be graphed. 
By default the 6 arrays described in `Description' are plotted. Useful
to produce individual postscript files of each array. For example,
specifying plotn = 1 in conjunction will postscript("01.ps") will
yield an array of plots of the quantile regression estimated coefficients.
}
\item{bcolor}{
color of the confidence band by default "gray".
}
\item{\dots}{
other optional arguments passed to `plot'.
}
}
\value{
Generates plots of object of class `khmal'. Please refer to
"Description" for further details.
}
\examples{
data(barro)
fit.Lonly _ rq.test.khmal(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, 
data = barro, location.scale = FALSE)
par(ask=interactive())
plot(fit.Lonly, var.list=c(2,4)) 
}
\keyword{htest}
\references{
Koenker, Roger and Zhijie Xiao (2000), "Inference on the Quantile
Regression Process'', unpublished.
\url{http://www.econ.uiuc.edu/~roger/research/inference/inference.html}
}
