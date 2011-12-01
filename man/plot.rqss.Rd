\name{plot.rqss}
\alias{plot.rqss}
\alias{plot.qss1}
\alias{plot.qss2}
\alias{plot.summary.rqss}
\title{Plot Method for rqss Objects}
\description{
Takes a fitted \code{rqss} object produced by \code{rqss()} and plots
the component smooth functions that make up the ANOVA decomposition.
Since the components "omit the intercept" the estimated intercept is added back
in -- this facilitates the comparison of quantile fits particularly.
For models with a partial linear component or several \code{qss} components
it may be preferable to plot the output of \code{predict.rqss}.
Note that these functions are intended to plot \code{rqss} objects only, attempting
to plot {summary.rqss}  objects just generates a warning message.
}
\usage{
\method{plot}{rqss}(x, rug = TRUE, jit = TRUE, bands = NULL, coverage = 0.95,
	add = FALSE, shade = TRUE, select = NULL, pages = 0, titles = NULL, \dots)
\method{plot}{qss1}(x, rug = TRUE, jit = TRUE, add = FALSE, \dots)
\method{plot}{qss2}(x, render = "contour", ncol = 100, zcol = NULL, \dots)
\method{plot}{summary.rqss}(x, ...)
}
\arguments{
  \item{x}{a fitted \code{rqss} object produced by \code{\link{rqss}()}.}
  \item{...}{additional arguments for the plotting algorithm}
  \item{rug}{if TRUE, a rugplot for the x-coordinate is plotted}
  \item{jit}{if TRUE, the x-values of the rug plot are jittered}
  \item{bands}{if TRUE, confidence bands for the smoothed effects are plotted, if
    "uniform" then uniform bands are plotted, if "both" then both the uniform
    and the pointwise bands are plotted.}
  \item{coverage}{desired coverage probability of confidence bands, if requested}
  \item{select}{vector of indices of qss objects to be plotted, by default all}
  \item{pages}{number of pages desired for the plots}
  \item{render}{a character specifying the rendering for bivariate fits;
    either \code{"contour"} (default) or \code{"rgl"}.  The latter
    requires package \pkg{rgl}.}
  \item{add}{if TRUE then add qss curve to existing (usually) scatterplot,
    otherwise initiate a new plot}
  \item{shade}{if TRUE then shade the confidence band}
  \item{titles}{title(s) as vector of character strings, by default titles are chosen for
	each plot as "Effect of CovariateName"}
  \item{ncol, zcol}{Only for \code{render = "rgl"}: number of colors and
    z values for color construction.}%% FIXME: be more specific
}
\value{
    The function produces plots for the ANOVA components as a side effect. For
\code{"qss1"} the \code{"add = TRUE"} can be used to overplot the fit on a 
scatterplot.  When there are multiple pages required \code{"par(ask = TRUE)"}
is turned on so that the plots may be examined sequentially.  If \code{bands != NULL}
then a list with three components for each qss component is returned (invisibly):
	\item{x}{The x coordinates of the confidence bands}
	\item{blo}{The y coordinates of the lower confidence curve, if 
		\code{bands = "both"} then this is a  matrix with two columns}
	\item{bhi}{The y coordinates of the upper confidence curve, if
	 	\code{bands = "both"} then this is a  matrix with two columns}	
}
\details{
If \code{bands == "uniform"} then the bands are uniform bands based on the 
Hotelling (1939) tube approach.  See also Naiman (1986), 
Johansen and Johnstone (1990), Sun and Loader (1994), 
and Krivobokova, Kneib, and Claeskens (2009), in particular the computation of
the "tube length" is based on the last of these references.  If \code{bands}
is non null, and not "uniform" then pointwise bands are returned.
Since bands for bivariate components are not (yet) supported, if requested
such components will be returned as \code{NULL}.
}
\author{ Roger Koenker }
\seealso{ \code{\link{rqss}}}
\references{
[1] Hotelling, H.  (1939): ``Tubes and Spheres in $n$-spaces, and a class
  of statistical problems,'' \emph{Am J. Math}, 61, 440--460.

[2] Johansen, S.,   and I.M. Johnstone  (1990): ``Hotelling's
  Theorem on the Volume of Tubes: Some Illustrations in Simultaneous
  Inference and Data Analysis,'' \emph{The Annals of Statistics}, 18, 652--684.

[3] Naiman, D. (1986)  Conservative confidence bands in curvilinear regression,  
  \emph{The Annals of Statistics}, 14, 896--906.

[4] Sun, J. and C.R. Loader, (1994) Simultaneous confidence bands for linear
regression and smoothing, \emph{The Annals of Statistics}, 22, 1328--1345.

[5] Krivobokova, T., T. Kneib, and G. Claeskens (2009) Simultaneous Confidence
Bands for Penalized Spline Estimators, preprint.

[6] Koenker, R. (2010) Additive Models for Quantile Regression:  Model Selection
and Confidence Bandaids, preprint.
}
\examples{
n <- 200
x <- sort(rchisq(n,4))
z <- x + rnorm(n)
y <- log(x)+ .1*(log(x))^2 + log(x)*rnorm(n)/4 + z
plot(x,y-z)
fN <- rqss(y~qss(x,constraint="N")+z)
plot(fN)
fI <- rqss(y~qss(x,constraint="I")+z)
plot(fI,  col="blue")
fCI <- rqss(y~qss(x,constraint="CI")+z)
plot(fCI, col="red")
%%keep objects for inspection : do not rm(n,x,y,z,fit)

## A bivariate example
data(CobarOre)
fCO <- rqss(z~qss(cbind(x,y),lambda=.08), data = CobarOre)
plot(fCO)
}
\keyword{regression}
\keyword{smooth}
\keyword{iplot}
