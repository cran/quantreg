\name{rq}
\alias{rq}
\title{
Quantile Regression 
}
\description{
Returns an object of class \code{"rq"} or \code{"rq.process"} that represents 
a quantile regression fit. 
}
\usage{
rq(formula, tau=.5, data, weights, na.action,
   method="br", contrasts, \dots) 
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a \code{~} operator, 
    and the terms, separated by \code{+} operators, on the right. 
  }
  \item{tau}{
    the quantile to be estimated, this is generally a number between 0 and 1, 
    but if specified outside this range, it is presumed that the solutions 
    for all values of \code{tau} in (0,1) are desired.  In the former case an
    object of class \code{"rq"} is returned, in the latter,
    an object of class \code{"rq.process"} is returned.
  }
  \item{data}{
    a data.frame in which to interpret the variables 
    named in the formula, or in the subset and the weights argument. 
    If this is missing, then the variables in the formula should be on the 
    search list.  This may also be a single number to handle some special  
    cases -- see below for details.   
  }
  \item{weights}{
    vector of observation weights; if supplied, the algorithm fits
    to minimize the sum of the weights multiplied into the
    absolute residuals. The length of weights must be the same as
    the number of observations.  The weights must be nonnegative
    and it is strongly recommended that they be strictly positive,
    since zero weights are ambiguous. 
  }
  \item{na.action}{
    a function to filter missing data. 
    This is applied to the model.frame after any subset argument has been used. 
    The default (with \code{na.fail}) is to create an error if any missing values are  
    found.  A possible alternative is \code{na.omit}, which 
    deletes observations that contain one or more missing values. 
  }
  \item{method}{
    the algorithmic method used to compute the fit.  There are currently 
    four options:   The default method is the modified  version of the
    Barrodale and Roberts algorithm for \eqn{l_1}{l1}-regression,
    used by \code{l1fit} in S, and is described in detail in 
    Koenker and d'Orey(1987, 1994),  default = \code{"br"}. 
    This is quite efficient for problems up to several thousand observations, 
    and may be used to compute the full quantile regression process.  It 
    also implements a scheme for computing confidence intervals for 
    the estimated parameters, based on inversion of a rank test described 
    in Koenker(1994).  For larger problems it is advantagous to use 
    the Frisch--Newton interior point method \code{"fn"}. 
    And very large problems one can use the Frisch--Newton approach after 
    preprocessing \code{"pfn"}.  Both of the latter methods are
    described in detail in Portnoy and Koenker(1997).   Finally, there
    is a fourth option \code{"fnc"} that enables the user to specify
    linear inequality constraints on the fitted coefficients; in this
    case one needs to specify the matrix \code{R} and the vector \code{r}
    representing the constraints in the form $Rb \geq r$.  See the
    examples 
  }
  \item{contrasts}{
    a list giving contrasts for some or all of the factors 
    default = \code{NULL} appearing in the model formula. 
    The elements of the list should have the same name as the variable 
    and should be either a contrast matrix (specifically, any full-rank 
    matrix with as many rows as there are levels in the factor), 
    or else a function to compute such a matrix given the number of levels. 
  }
  \item{...}{
    additional arguments for the fitting routines 
    (see \code{\link{rq.fit.br}} and \code{\link{rq.fit.fn}}
    and the functions they call). 
  }
}
\value{
  See \code{\link{rq.object}} and \code{\link{rq.process.object}} for details. 
}
\examples{
data(stackloss)
rq(stack.loss ~ stack.x,.5)  #median (l1) regression  fit for the stackloss data. 
rq(stack.loss ~ stack.x,.25)  #the 1st quartile, 
        #note that 8 of the 21 points lie exactly on this plane in 4-space! 
rq(stack.loss ~ stack.x, tau=-1)   #this returns the full rq process
rq(rnorm(50) ~ 1, ci=FALSE)    #ordinary sample median --no rank inversion ci
rq(rnorm(50) ~ 1, weights=runif(50),ci=FALSE)  #weighted sample median 
#plot of engel data and some rq lines see KB(1982) for references to data
data(engel)
attach(engel)
plot(x,y,xlab="household income",ylab="food expenditure",cex=.5)
taus <- c(.05,.1,.25,.75,.9,.95)
xx <- seq(min(x),max(x),100)
for(tau in taus){
        f <- rq((y)~(x),tau=tau)$coef[,1]
        yy <- (f[1]+f[2]*(xx))
        lines(xx,yy)
        }
#Example to illustrate inequality constrained fitting
n <- 100
p <- 5
X <- matrix(rnorm(n*p),n,p)
y <- .95*apply(X,1,sum)+rnorm(n)
#constrain slope coefficients to lie between zero and one
R <- cbind(0,rbind(diag(p),-diag(p)))
r <- c(rep(0,p),-rep(1,p))
rq(y~X,R=R,r=r,method="fnc")
}
\section{Method}{
The function computes an estimate on the tau-th conditional quantile
function of the response, given the covariates, as specified by the
formula argument.  Like \code{lm()}, the function presumes a linear
pecification for the quantile regression model, i.e. that the formula
defines a model that is linear in parameters.  For non-linear quantile
regression see the package \code{nlrq()}.  
The function minimizes a weighted sum of absolute
residuals that can be formulated as a linear programming problem.  As
noted above, there are three different algorithms that can be chosen
depending on problem size and other characteristics.  For moderate sized
problems (\eqn{n \ll 5,000, p \ll 20}{n << 5,000, p << 20}) it is recommended that the default
\code{"br"} method be used. There are several choices of methods for
computing confidence intervals and associated test statistics.  Using
\code{"br"} the default approach produces confidence intervals for each
of the estimated model parameters based on inversion of a rank test.
See the documentation for \code{\link{rq.fit.br}} for further details
and options.  For larger problems, the \code{"fn"} and \code{"pfn"} are
preferred, and there are several methods of computing standard errors
and associated test statistics described in the help files for
\code{\link{rq.fit.fn}}, and \code{\link{summary.rq}}.
}

\keyword{regression}
\references{
[1] Koenker, R. W. and Bassett, G. W. (1978). Regression quantiles, 
\emph{Econometrica}, \bold{46}, 33--50. 

[2] Koenker, R.W. and d'Orey (1987, 1994). Computing regression quantiles. 
\emph{Applied Statistics}, \bold{36}, 383--393, and \bold{43}, 410--414. 

[3] Gutenbrunner, C. Jureckova, J. (1991). 
Regression quantile and regression rank score process in the 
linear model and derived statistics, \emph{Annals of Statistics},
\bold{20}, 305--330.

[4] Koenker, R. W. (1994). Confidence Intervals for regression quantiles, in 
P. Mandl and M. Huskova (eds.), \emph{Asymptotic Statistics}, 349--359,  
Springer-Verlag, New York.   

[5] Koenker, R. and S. Portnoy (1997) The Gaussian Hare and the Laplacean 
Tortoise:  Computability of Squared-error vs Absolute Error Estimators, 
(with discussion).  \emph{Statistical Science,} \bold{12}, 279-300.

There is also recent information available at the URL:
\url{http://www.econ.uiuc.edu}.
}
\seealso{
  \code{\link{summary.rq}}, \code{\link{rq.object}},
  \code{\link{rq.process.object}}
}
% Converted by Sd2Rd version 0.3-3.
