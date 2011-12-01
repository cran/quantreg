\name{predict.rqss}
\alias{predict.rqss}
\alias{predict.qss1}
\alias{predict.qss2}
\title{Predict from fitted nonparametric quantile regression smoothing spline models}
\description{
  Additive models for nonparametric quantile regression using total
  variation penalty methods can be fit with the \code{\link{rqss}}
  function.  Univarariate and bivariate components can be predicted
  using these functions.
}
\usage{
\method{predict}{rqss}(object, newdata, interval = "none", level = 0.95, ...)
\method{predict}{qss1}(object, newdata, ...)
\method{predict}{qss2}(object, newdata, ...)
}
\arguments{
\item{object}{ is a fitted object produced by \code{\link{rqss}} }
\item{newdata}{ a data frame describing the observations at which
  prediction is to be made.  For qss components, newdata should
  lie in strictly within the convex hull of the fitting data.  Newdata
  corresponding to the partially linear component of the model
  may require caution concerning the treatment of factor levels, if any.}
\item{interval}{If set to \code{confidence} then a \code{level} confidence interval
for the predictions is returned.}
\item{level}{intended coverage probability for the confidence intervals}
\item{\dots}{ optional arguments }
}
\details{
For both univariate and bivariate prediction linear interpolation is
done.  In the bivariate case, this involves computing barycentric
coordinates of the new points relative to their enclosing triangles.
It may be of interest to plot individual components of fitted rqss
models:  this is usually best done by fixing the values of other
covariates at reference values typical of the sample data and
predicting the response at varying values of one qss term at a
time.   Direct use of the \code{predict.qss1} and \code{predict.qss2} functions
is discouraged since it usually corresponds to predicted values
at absurd  reference values of the other covariates, i.e. zero.
}
\value{
A vector of predictions, or in the case that \code{interval = "confidence")}
a matrix whose first column is the vector of predictions and whose second and
third columns are the lower and upper confidence limits for each prediction.
}
\author{ R. Koenker }
\seealso{ \code{\link{rqss}}
}
\examples{
n <- 200
lam <- 2
x <- sort(rchisq(n,4))
z <- exp(rnorm(n)) + x
y <- log(x)+ .1*(log(x))^2 + z/4 +  log(x)*rnorm(n)/4
plot(x,y - z/4 + mean(z)/4)
Ifit <- rqss(y ~ qss(x,constraint="I") + z)
sfit <- rqss(y ~ qss(x,lambda = lam) + z)
xz <- data.frame(z = mean(z),
                 x = seq(min(x)+.01,max(x)-.01,by=.25))
lines(xz[["x"]], predict(Ifit, xz), col=2)
lines(xz[["x"]], predict(sfit, xz), col=3)
legend(10,2,c("Increasing","Smooth"),lty = 1, col = c(2,3))
title("Predicted Median Response at Mean Value of z")
%%keep objects for inspection : do not rm(x,y,z,xz,fit)

## Bivariate example -- loads pkg "tripack"
data(CobarOre)
fit <- rqss(z ~ qss(cbind(x,y), lambda=.08),
            data= CobarOre)
plot(fit, col="grey",
     main = "CobarOre data -- rqss(z ~ qss(cbind(x,y)))")
T <- with(CobarOre, tri.mesh(x, y))
set.seed(77)
ndum <- 100
xd <- with(CobarOre, runif(ndum, min(x), max(x)))
yd <- with(CobarOre, runif(ndum, min(y), max(y)))
table(s <- in.convex.hull(T, xd, yd))
pred <- predict(fit, data.frame(x = xd[s], y = yd[s]))
contour(interp(xd[s],yd[s], pred),
        col="red", add = TRUE)
}
\keyword{regression}
\keyword{smooth}
\keyword{robust}
