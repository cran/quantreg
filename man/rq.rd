\name{rq}
\alias{rq}
\alias{rq.formula}
\title{
Quantile Regression
}
\description{
Perform a quantile regression on a design matrix, x, of explanatory variables and a vector, y, of responses. 
}
\usage{
rq(x, y, tau=-1, alpha=.1, dual=TRUE, int=TRUE, tol=1e-4, ci = TRUE, 
         method="score", interpolate=TRUE, tcrit=TRUE, hs=TRUE)
rq.formula(formula, data=list(), subset, na.action, tau=-1, 
         alpha = 0.10000000000000001, dual = TRUE, 
         tol = 0.0001, ci = TRUE, method="score", interpolate = TRUE, 
         tcrit = TRUE, hs=TRUE)

}
\arguments{
\item{x}{
vector or matrix of explanatory variables.  If  a  matrix,
each  column represents a variable and each row represents
an observation (or case).  This should not contain  column
of  1s unless the argument intercept is FALSE.  The number
of rows of x should equal the number of elements of  y,  and
there  should  be fewer columns than rows.  
If x is missing, rq() computes the ordinary
sample quantile(s) of y.
}
\item{y}{
response vector with as many observations as the number of rows of x. 
}
\item{tau}{
desired quantile. If tau is missing or outside the range [0,1]
then all the regression quantiles are computed and the corresponding primal and dual solutions are returned.
}
\item{alpha}{
level of significance for the confidence intervals; default is set at 10\%.
}
\item{dual}{
return the dual solution if TRUE (default).
}
\item{int}{
flag for intercept; if TRUE (default) an intercept term is included in the regression.
}
\item{tol}{
tolerance parameter for rq computations.
}
\item{ci}{
flag for confidence interval; if TRUE (default) the confidence intervals are 
returned.
}
\item{method}{
if method="score" (default), ci is computed using regression rank score inversion;
if method="sparsity", ci is computed using sparsity function.
}
\item{interpolate}{
if TRUE (default), the smoothed confidence intervals are returned.
}
\item{tcrit}{
if tcrit=T (default), a finite sample adjustment of the critical point is
performed using Student's t quantile, else the standard Gaussian quantile is 
used.
}
\item{hs}{
logical flag to use Hall-Sheather's sparsity estimator (default); otherwise Bofinger's
version is used.
}}
\value{


\item{coef}{
the estimated parameters of the tau-th conditional quantile function.
}
\item{resid}{
the estimated residuals of the tau-th conditional quantile function.
}
\item{dual}{
the dual solution (if dual=T).
}
\item{h}{
the index of observations in the basis.
}
\item{ci}{
confidence intervals (if ci=T).
}}
\value{


\item{sol}{
a  (p+2) by m matrix whose first row contains the 'breakpoints'
tau_1,tau_2,\dots{}tau_m, of the quantile function, 
i.e. the values in [0,1] at which the
solution changes, row two contains the corresponding quantiles
evaluated at the mean design point, i.e. the inner product of
xbar and b(tau_i), and the last p rows of the matrix give b(tau_i).
The solution b(tau_i) prevails from tau_i to tau_i+1.
}
\item{dsol}{
the matrix of dual solutions corresponding to the primal solutions in sol.
This is an n by m matrix whose ij-th entry is 1 if y_i > x_i b(tau_j), 
is 0 if y_i < x_i b(tau_j),  and is between 0 and 1 otherwise, i.e. if
the residual is zero.  See Gutenbrunner and Jureckova(1991) for a
detailed discussion of the statistical interpretation of dsol.
}
\item{h}{
the matrix of observations indices in the basis corresponding to sol or dsol.
}}
\author{Roger Koenker, \email{roger@ysidro.econ.uiuc.edu}, 
        \url{http://www.econ.uiuc.edu/~roger/research/rq/rq.html}. Ported
        to R, and added rq.formula, by Kjetil Halvorsen.}
\examples{
data(stackloss)
rq(stack.x, stack.loss, .5)  #the l1 estimate for the stackloss data
rq(stack.x, stack.loss, tau=.5, ci=T, method="score")  #same as above with 
	#regression rank score inversion confidence interval
rq(stack.x, stack.loss, .25)  #the 1st quartile, 
	#note that 8 of the 21 points lie exactly 
	#on this plane in 4-space
rq(stack.x, stack.loss, -1)   #this gives all of the rq solutions
rq(y=rnorm(10), method="sparsity")	#ordinary sample quantiles
\dontrun{data(Patacamaya)}               # an example with formula
\dontrun{ z0.1 <- rq.formula(y ~ a+tipo, data=Patacamaya, na.action=na.omit, tau=0.1)}
\dontrun{z0.1$coef}
\dontrun{z0.1$ci}
}
\section{METHOD}{
The algorithm used is a modification of the Barrodale and Roberts
algorithm for l1-regression, l1fit in S, and is described in detail
in Koenker and d"Orey(1987).
}
\keyword{regression}
\references{
[1] Koenker, R.W. and Bassett, G.W. (1978). Regression quantiles, Econometrica, 46, 33-50.



[2] Koenker, R.W. and d'Orey (1987). Computing Regression Quantiles. Applied Statistics, 36, 383-393.



[3] Gutenbrunner, C. Jureckova, J. (1991). 
Regression quantile and regression rank score process in the 
linear model and derived statistics, Annals of Statistics, 20, 305-330.



[4] Koenker, R.W. and d'Orey (1994).  Remark on Alg. AS 229: Computing Dual
Regression Quantiles and Regression Rank Scores, Applied Statistics, 43, 410-414.



[5] Koenker, R.W. (1994). Confidence Intervals for Regression Quantiles, in
P. Mandl and M. Huskova (eds.), Asymptotic Statistics, 349-359, Springer-Verlag,
New York.



}
\section{SEE ALSO}{
trq and qrq for further details and references.


}
% Converted by Sd2Rd version 0.3-2.
