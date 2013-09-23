\name{dynrq}
\alias{dynrq}
\alias{print.dynrq}
\alias{print.dynrqs}
\alias{summary.dynrq}
\alias{summary.dynrqs}
\alias{print.summary.dynrq}
\alias{print.summary.dynrqs}
\alias{time.dynrq}
\alias{index.dynrq}
\alias{start.dynrq}
\alias{end.dynrq}

\title{Dynamic Linear Quantile Regression}

\description{
Interface to \code{\link{rq.fit}} and \code{\link{rq.wfit}} for fitting dynamic linear 
quantile regression models.  The interface is based very closely
on Achim Zeileis's \code{dynlm} package.  In effect, this is  mainly
``syntactic sugar'' for formula processing, but one should never underestimate
the value of good, natural sweeteners.
}

\usage{dynrq(formula, tau = 0.5, data, subset, weights, na.action, method = "br",
  contrasts = NULL, start = NULL, end = NULL, ...)}

\arguments{
  \item{formula}{a \code{"formula"} describing the linear model to be fit.
    For details see below and \code{\link{rq}}.}
  \item{tau}{the quantile(s) to be estimated, may be vector valued, but all
    all values must be in (0,1).} 
  \item{data}{an optional \code{"data.frame"} or time series object (e.g.,
    \code{"ts"} or \code{"zoo"}), containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from which
    \code{rq} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{weights}{an optional vector of weights to be used
    in the fitting process. If specified, weighted least squares is used
    with weights \code{weights} (that is, minimizing \code{sum(w*e^2)});
    otherwise ordinary least squares is used.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}. Another possible value is
    \code{NULL}, no action. Note, that for time series regression
    special methods like \code{\link{na.contiguous}}, \code{\link[zoo]{na.locf}}
    and \code{\link[zoo]{na.approx}} are available.}
  \item{method}{the method to be used; for fitting, by default
    \code{method = "br"} is used; \code{method = "fn"} employs
    the interior point (Frisch-Newton) algorithm.  The latter is advantageous
    for problems with sample sizes larger than about 5,000.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
  \item{start}{start of the time period which should be used for fitting the model.}
  \item{end}{end of the time period which should be used for fitting the model.}
  \item{\dots}{additional arguments to be passed to the low level
    regression fitting functions.}
}

\details{
The interface and internals of \code{dynrq} are very similar to \code{\link{rq}},
but currently \code{dynrq} offers two advantages over the direct use of
\code{rq} for time series applications of quantile regression: 
extended formula processing, and preservation of time series attributes.  
Both features have been shamelessly lifted from Achim Zeileis's
package \code{"dynlm"}.

For specifying the \code{formula} of the model to be fitted, there are several
functions available which allow for convenient specification
of dynamics (via \code{d()} and \code{L()}) or linear/cyclical patterns
(via \code{trend()}, \code{season()}, and \code{harmon()}).
These new formula functions require that their arguments are time
series objects (i.e., \code{"ts"} or \code{"zoo"}).

Dynamic models: An example would be \code{d(y) ~ L(y, 2)}, where
\code{d(x, k)} is \code{diff(x, lag = k)} and \code{L(x, k)} is
\code{lag(x, lag = -k)}, note the difference in sign. The default
for \code{k} is in both cases \code{1}. For \code{L()}, it
can also be vector-valued, e.g., \code{y ~ L(y, 1:4)}. 

Trends: \code{y ~ trend(y)} specifies a linear time trend where
\code{(1:n)/freq} is used by default as the covariate, \code{n} is the 
number of observations and \code{freq} is the frequency of the series
(if any, otherwise \code{freq = 1}). Alternatively, \code{trend(y, scale = FALSE)}
would employ \code{1:n} and \code{time(y)} would employ the original time index.

Seasonal/cyclical patterns: Seasonal patterns can be specified
via \code{season(x, ref = NULL)} and harmonic patterns via
\code{harmon(x, order = 1)}.  \code{season(x, ref = NULL)} creates a factor 
with levels for each cycle of the season. Using
the \code{ref} argument, the reference level can be changed from the default
first level to any other. \code{harmon(x, order = 1)} creates a matrix of
regressors corresponding to \code{cos(2 * o * pi * time(x))} and 
\code{sin(2 * o * pi * time(x))} where \code{o} is chosen from \code{1:order}.

See below for examples. 

Another aim of \code{dynrq} is to preserve 
time series properties of the data. Explicit support is currently available 
for \code{"ts"} and \code{"zoo"} series. Internally, the data is kept as a \code{"zoo"}
series and coerced back to \code{"ts"} if the original dependent variable was of
that class (and no internal \code{NA}s were created by the \code{na.action}).

}

\seealso{\code{\link[zoo]{zoo}}, 
\code{\link[dynlm]{dynlm}},
\code{\link[zoo]{merge.zoo}}}

\examples{
###########################
## Dynamic Linear Quantile Regression Models ##
###########################

## multiplicative median SARIMA(1,0,0)(1,0,0)_12 model fitted to UK seatbelt data
     data("UKDriverDeaths", package = "datasets")
     uk <- log10(UKDriverDeaths)
     dfm <- dynrq(uk ~ L(uk, 1) + L(uk, 12))
     dfm

     dfm3 <- dynrq(uk ~ L(uk, 1) + L(uk, 12),tau = 1:3/4)
     summary(dfm3)
 ## explicitly set start and end
     dfm1 <- dynrq(uk ~ L(uk, 1) + L(uk, 12), start = c(1975, 1), end = c(1982, 12))
 ## remove lag 12
     dfm0 <- update(dfm1, . ~ . - L(uk, 12))
     tuk1  <- anova(dfm0, dfm1)
 ## add seasonal term
     dfm1 <- dynrq(uk ~ 1, start = c(1975, 1), end = c(1982, 12))
     dfm2 <- dynrq(uk ~ season(uk), start = c(1975, 1), end = c(1982, 12))
     tuk2 <- anova(dfm1, dfm2)
 ## regression on multiple lags in a single L() call
     dfm3 <- dynrq(uk ~ L(uk, c(1, 11, 12)), start = c(1975, 1), end = c(1982, 12))
     anova(dfm1, dfm3)

###############################
## Time Series Decomposition ##
###############################

## airline data
data("AirPassengers", package = "datasets")
ap <- log(AirPassengers)
fm <- dynrq(ap ~ trend(ap) + season(ap), tau = 1:4/5)
sfm <- summary(fm)
plot(sfm)

## Alternative time trend specifications:
##   time(ap)                  1949 + (0, 1, ..., 143)/12
##   trend(ap)                 (1, 2, ..., 144)/12
##   trend(ap, scale = FALSE)  (1, 2, ..., 144)

}

\keyword{regression}
