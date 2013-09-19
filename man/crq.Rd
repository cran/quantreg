\name{crq}
\alias{crq}
\alias{crq.fit.por}
\alias{crq.fit.pow}
\alias{crq.fit.pen}
\alias{print.crq}
\alias{print.crq}
\alias{coef.crq}
\alias{predict.crq}
\alias{predict.crqs}
\alias{Curv}
\title{Functions to fit censored quantile regression models} 
\description{Fits a conditional quantile regression model for censored data. There
are three distinct methods:  the first is the fixed censoring method
of Powell (1986) as implemented by Fitzenberger (1996), the second is the random
censoring method of Portnoy (2003).  The third method is based on Peng and Huang (2008).} 
\usage{
crq(formula, taus, data, subset, weights, na.action, 
	method = c("Powell", "Portnoy", "PengHuang"), contrasts = NULL, ...)
crq.fit.pow(x, y, yc, tau=0.5, weights=NULL, start, left=TRUE, maxit = 500)
crq.fit.pen(x, y, cen, weights=NULL, grid, ctype = "right") 
crq.fit.por(x, y, cen, weights=NULL, grid, ctype = "right") 
Curv(y, yc, ctype=c("left","right"))
\method{print}{crq}(x, ...)
\method{print}{crq}(x, ...)
\method{predict}{crq}(object, newdata,  ...)
\method{predict}{crqs}(object, newdata, type = NULL, ...)
\method{coef}{crq}(object,taus = 1:4/5,...)
}
\arguments{
  \item{formula}{A formula object, with the response on the left of the `~'
        operator, and the terms on the right.  The response must be a
        \code{Surv} object as returned by either the \code{Curv} or \code{Surv} 
	function. For the Powell method, the Surv object should
	be created by \code{Curv} and have arguments (event time, censoring time,type), 
	where "type" can take values either "left" or "right". 
	The default (for historical reasons) for type in this case is "left".
	For the Portnoy and Peng and Huang  methods the \code{Surv}  should be created 
	with the usual \code{Surv} function and have (event time, censoring indicator).}  
  \item{y}{The event time.}
  \item{newdata}{An optional data frame in which to look for variables with which 
	to predict. If omitted, the fitted values are used.}
  \item{grid}{A vector of taus on which the quantile process should be evaluated. This
	should be monotonic, and take values in (0,1).  For the "Portnoy"  method,
	grid = "pivot" computes the full solution for all distinct taus.  The "Portnoy"
	method also enforces an equally spaced grid, see the code for details.}
  \item{x}{An object of class \code{crq} or \code{crq}.}
  \item{object}{An object of class \code{crq} or \code{crq}.}
  \item{yc}{The censoring times for the "Powell" method.}
  \item{ctype}{Censoring type: for the "Powell" method, used in \code{Curv}, by 
	default "left".  If you don't like "left", maybe you will like "right".
	Note that for fixed censoring assumed in the "Powell" method, censoring
	times \code{yc} must be provided for all observations and the event
	times \code{y} must satisfy the (respective) inequality constraints.
	For the Portnoy and Peng-Huang methods ctype is determined by the
	specification of the response as specified in \code{Surv}.
        }
  \item{type}{specifies either "left" or "right" as the form of censoring 
	 in the \code{Surv} function for the "Portnoy" and "PengHuang"  methods.}
  \item{cen}{The censoring indicator for the "Portnoy" and "PengHuang"  methods.}
  \item{maxit}{Maximum number of iterations allowed for the "Powell" methods.}
  \item{start}{The starting value for the coefs for the "Powell" method.  Because
	the Fitzenberger algorithm stops when it achieves a local minimum
	of the Powell objective function, the starting value acts as an
	a priori "preferred point".  This is advantageous in some instances
	since the global Powell solution can be quite extreme. By default the
	starting value is the "naive rq" solution that  treats all the censored
	observations as uncensored.  If \code{start} is equal to "global"
	then an attempt is made to compute to global optimum of the Powell
	objective.  This entails an exhaustive evaluation of all n choose p
	distinct basic solution so is rather impractical for moderately large
	problems. Otherwise, the starting value can specify a set of p indices
	from 1:n defining an initial basic solution, or it may specify a p-vector
	of initial regression coefficients.  In the latter case the initial basic
	solution is the one closest to the specified parameter vector.}
  \item{left}{A logical indicator for left censoring for the "Powell" method.}
  \item{taus}{The quantile(s) at which the model is to be estimated.}
  \item{tau}{The quantile at which the model is to be estimated.}
  \item{data}{A data.frame in which to interpret the variables named in the
          `formula',  in the `subset', and the `weights' argument.}
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.}
  \item{weights}{vector of observation weights; if supplied, the algorithm
          fits to minimize the sum of the weights multiplied into the
          absolute residuals. The length of weights vector must be the same as
          the number of observations.  The weights must be nonnegative
          and it is strongly recommended that they be strictly
          positive, since zero weights are ambiguous.} 
  \item{na.action}{a function to filter missing data.  This is applied to the
          model.frame after any subset argument has been used.  The
          default (with 'na.fail') is to create an error if any missing
          values are   found.  A possible alternative is 'na.omit',
          which  deletes observations that contain one or more missing
          values. } 
  \item{method}{The method used for fitting.  There are currently
	two options: method "Powell" computes the Powell estimator using
	the algorithm of Fitzenberger (1996), method "Portnoy" computes the
	Portnoy (2003) estimator.  The  method is "PengHuang" uses the method
	of Peng and Huang (2007), in this case the variable "grid"
	can be passed to specify the vector of quantiles at which the solution
	is desired.}
  \item{contrasts}{a list giving contrasts for some or all of the factors 
          default = 'NULL' appearing in the model formula.  The
          elements of the list should have the same name as the
          variable  and should be either a contrast matrix
          (specifically, any full-rank  matrix with as many rows as
          there are levels in the factor),  or else a function to
          compute such a matrix given the number of levels.} 
  \item{...}{additional arguments for the fitting routine, for method "Powell"
	it may be useful to pass starting values of the regression parameter
	via the argument "start", while for methods "Portnoy" or "PengHuang"
	one may wish to specify an alternative to the default grid for evaluating 
	the fit.}
}
\details{The Fitzenberger algorithm uses a variant of the Barrodale and Roberts
	simplex method.  Exploiting the fact that the solution must be characterized
	by an exact fit to p points when there are p parameters to be estimated,
	at any trial basic solution it computes the directional derivatives in the 
	2p distinct directions
	and choses the direction that (locally) gives steepest descent.  It then
	performs a one-dimensional line search to choose the new basic observation
	and continues until it reaches a local mimumum.  By default it starts at
	the naive \code{rq} solution ignoring the censoring;  this has the (slight)
	advantage that the estimator is consequently equivariant to canonical
	transformations of the data.  Since the objective function is no longer convex
	there can be no guarantee that this  produces a global minimum estimate.
	In small problems exhaustive search over solutions defined by p-element
	subsets of the n observations can be used, but this quickly becomes
	impractical for large p and n.   This global version of the Powell
	estimator can be invoked by specifying \code{start = "global"}. Users
	interested in this option would be well advised to compute \code{choose(n,p)}
	for their problems before trying it.  The method operates by pivoting
	through this many distinct solutions and choosing the one that gives the
	minimal Powell objective.  The algorithm used for the Portnoy 
	method is described in considerable detail in Portnoy (2003). 
	Both the Portnoy and Peng-Huang estimators may be unable to compute
	estimates of the conditional quantile parameters in the upper tail of
	distribution.  Like the Kaplan-Meier estimator, when censoring is heavy
	in the upper tail the estimated distribution is defective  and quantiles
	are only estimable on a sub-interval of (0,1).
	The Peng and Huang estimator can be
	viewed as a generalization of the Nelson Aalen estimator of the cumulative
	hazard function,  and can be formulated as a variant of the conventional
	quantile regression dual problem.  See Koenker (2008) for further details.
	This paper is available from the package with \code{vignette("crq")}.}
\value{An object of class \code{crq}.}

\references{

Fitzenberger, B.  (1996): ``A Guide to Censored Quantile
Regressions,'' in \emph{Handbook of Statistics}, ed. by C.~Rao,   and
G.~Maddala. North-Holland: New York.

Fitzenberger, B.  and P. Winker (2007): ``Improving the Computation of
Censored Quantile Regression Estimators,'' CSDA, 52,  88-108.

Koenker, R. (2008): ``Censored Quantile Regression Redux,'' \emph{J. 
Statistical Software}, 27, \url{http://www.jstatsoft.org/v27/i06}.

Peng, L and Y Huang, (2008) Survival Analysis with Quantile Regression Models,
\emph{J. Am. Stat. Assoc.}, 103, 637-649. 

Portnoy, S. (2003) ``Censored Quantile Regression,'' \emph{JASA},
98,1001-1012.

Powell, J. (1986) ``Censored Regression Quantiles,'' \emph{J.
Econometrics}, 32, 143--155.
}


\author{Steve Portnoy and  Roger Koenker}

\seealso{\code{\link{summary.crq}}}

\examples{
# An artificial Powell example
set.seed(2345)
x <- sqrt(rnorm(100)^2)
y <-  -0.5 + x +(.25 + .25*x)*rnorm(100)
plot(x,y, type="n")
s <- (y > 0)
points(x[s],y[s],cex=.9,pch=16)
points(x[!s],y[!s],cex=.9,pch=1)
yLatent <- y
y <- pmax(0,y)
yc <- rep(0,100)
for(tau in (1:4)/5){
        f <- crq(Curv(y,yc) ~ x, tau = tau, method = "Pow")
        xs <- sort(x)
        lines(xs,pmax(0,cbind(1,xs)\%*\%f$coef),col="red")
        abline(rq(y ~ x, tau = tau), col="blue")
        abline(rq(yLatent ~ x, tau = tau), col="green")
        }
legend(.15,2.5,c("Naive QR","Censored QR","Omniscient QR"),
        lty=rep(1,3),col=c("blue","red","green"))

# crq example with left censoring
set.seed(1968)
n <- 200
x <-rnorm(n)
y <- 5 + x + rnorm(n)
plot(x,y,cex = .5)
c <- 4 + x + rnorm(n)
d <- (y > c)
points(x[!d],y[!d],cex = .5, col = 2)
f <- crq(Surv(y,d,type = "left") ~ x, method = "Portnoy")
g <- summary(f)
for(i in 1:4) abline(coef(g[[i]])[,1])
}

\keyword{survival}
\keyword{regression}
