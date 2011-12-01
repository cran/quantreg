\name{predict.rq}
\alias{predict.rq}
\alias{predict.rqs}
\alias{predict.rq.process}
\title{Quantile Regression Prediction}
\description{
Prediction based on fitted quantile regression model
}
\usage{
\method{predict}{rq}(object, newdata, type = "none", interval = c("none", "confidence"), 
	level = .95, na.action = na.pass, ...)
\method{predict}{rqs}(object, newdata, type = "Qhat", stepfun = FALSE, na.action = na.pass, ...)
\method{predict}{rq.process}(object, newdata, type = "Qhat", stepfun = FALSE, na.action = na.pass, ...)
}
\arguments{
  \item{object}{ object of class rq or rqs or rq.process produced by \code{rq} }
  \item{newdata}{An optional data frame in which to look for variables with
    which to predict.  If omitted, the fitted values are used.}
  \item{interval}{type of interval desired:  default is  'none', when set to
    'confidence' the function returns a matrix predictions with point predictions
    for each of the 'newdata' points as well as lower and upper confidence limits.}
  \item{level}{converage probability for the 'confidence' intervals.}
  \item{type}{For \code{predict.rq}, the  method for 'confidence' intervals, if desired. 
    If 'percentile' then one of the bootstrap methods is used to generate percentile 
    intervals for each prediction, if 'direct' then a version of the Portnoy and Zhou 
    (1998) method is used, and otherwise an estimated covariance matrix for the parameter
    estimates is used.  Further arguments to determine the choice of bootstrap
    method or covariance matrix estimate can be passed via the \dots argument.
    For \code{predict.rqs} and \code{predict.rq.process} when \code{stepfun = TRUE},
    \code{type} is either "Qhat" or "Fhat" depending  on whether the user would
    like to have estimates of the conditional quantile or distribution functions
    respectively.  As noted below such estimates can be monotonized with the function
    \code{rearrange}.}
  \item{stepfun}{If 'TRUE' return stepfunctions otherwise return matrix of predictions.
    these functions can be estimates of either the conditional quantile or distribution
    functions depending upon the \code{type} argument.  When \code{stepfun = FALSE}
    a matrix of point estimates of the conditional quantile function at the points
    specified by the \code{newdata} argument. }
  \item{na.action}{ function determining what should be done with missing values
    in 'newdata'.  The default is to predict 'NA'.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
     Produces predicted values, obtained by evaluating the quantile
     regression function in the frame 'newdata' (which defaults to
     'model.frame(object)'.  These predictions purport to estimate
    the conditional quantile function of the response variable of
    the fitted model evaluated at the covariate values specified
    in "newdata" and the quantile(s) specified by the "tau" argument.
    Several methods are provided to compute confidence intervals for
    these predictions.  
}
\value{
     A vector or matrix of predictions, depending upon the setting of 
     'interval'.  In the case that there are multiple taus in \code{object}
     when object is of class 'rqs' setting 'stepfun = TRUE'  will produce a
     \code{stepfun} object or a list of \code{stepfun} objects.
     The function \code{rearrange} can be used to monotonize these
     step-functions, if desired.
}
\references{ 
Zhou, Kenneth Q. and Portnoy, Stephen L. (1998) 
Statistical inference on heteroscedastic models based on regression quantiles 
Journal of Nonparametric Statistics, 9, 239-260 
}
\author{R. Koenker}
\seealso{ \code{\link{rq}} \code{\link{rearrange}}}
\examples{
data(airquality)
airq <- airquality[143:145,]
f <- rq(Ozone ~ ., data=airquality)
predict(f,newdata=airq)
f <- rq(Ozone ~ ., data=airquality, tau=1:19/20)
fp <- predict(f, newdata=airq, stepfun = TRUE)
fpr <- rearrange(fp)
plot(fp[[2]],main = "Conditional Ozone Quantile Prediction")
lines(fpr[[2]], col="red")
legend(.2,20,c("raw","cooked"),lty = c(1,1),col=c("black","red"))
fp <- predict(f, newdata=airq, type = "Fhat", stepfun = TRUE)
fpr <- rearrange(fp)
plot(fp[[2]],main = "Conditional Ozone Distribution Prediction")
lines(fpr[[2]], col="red")
legend(20,.4,c("raw","cooked"),lty = c(1,1),col=c("black","red"))
  }
\keyword{regression}
