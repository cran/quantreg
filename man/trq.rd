\name{trq}
\alias{trq}
\title{
Function to compute analogues of the trimmed mean 
for the linear regression model.
}
\description{
The function returns 
a regression trimmed mean and some associated test statistics.
The proportion a1 is trimmed from the lower tail and a2 from the upper tail.
If a1+a2=1 then a result is returned for the a1 quantile.
If a1+a2<1 two methods of trimming are possible described below as
"primal" and "dual". The function "trq.print" may be used to print results in the style of ls.print.
}
\usage{
trq(x, y, a1=0.1, a2,  int=TRUE, z,  method="primal", tol=1e-4)
}
\arguments{
\item{x}{
vector or matrix of explanatory variables.  If  a  matrix,
each  column represents a variable and each row represents
an observation (or case).  This should not contain  column
of  1s unless the argument intercept is FALSE.  The number
of rows of x should equal the number of elements of  y,  and
there  should  be fewer columns than rows.  Missing values
are not  allowed.
}
\item{y}{
reponse vector with as many observations as the number of rows of x.  
Missing value are not allowed.
}
\item{a1}{
the lower trimming proportion; defaults to .1 if missing.
}
\item{a2}{
the upper trimming proportion; defaults to a1 if missing.
}
\item{int}{
flag for intercept; if TRUE, an intercept term is included in regression model. 
The default includes an intercept term.
}
\item{z}{
structure returned by the function 'rq' with tau <0 or >1. If missing, the 
function rq(x,y,int=int) is automatically called to generate this argument.
If several calls to trq are anticipated for the same data this avoids
recomputing the rq solution for each call.
}
\item{method}{
method to be used for the trimming. 
If the choice is "primal", as is the default, a trimmed mean 
of the primal regression quantiles  is computed based on the sol array 
in the 'rq' structure. 
If the method is "dual", a weighted least-squares fit is done 
using the dual solution in the 'rq' structure to construct weights.
The former method is discussed in detail in Koenker and Potnoy(1987)
the latter in Ruppert and Carroll(1980) and Gutenbrunner and Jureckova(1991).
}
\item{tol}{
Tolerance parameter for rq computions
}}
\value{


\item{coef}{
estimated coeficient vector
}
\item{resid}{
residuals from the fit.
}
\item{cov}{
the estimated covariance matrix for the coeficient vector. 
}
\item{v}{
the scaling factor of the covariance matrix under iid error assumption:
cov=v*(x'x)^(-1).
}
\item{wt}{
the weights used in the least squares computation,  Returned only when
method="dual".
}
\item{d}{
the bandwidth used to compute the sparsity function. 
Returned only when a1+a2=1.
}}
\examples{
x <- -10:10; y <- 0.2 * x + rt(x, df=3)
z <- rq(x,y)        #z gets the full regression quantile structure
trq(x,y, .05, z=z)  #5\% symmetric primal trimming # Error, which also occurs in S-Plus. 
trq(x,y, .01, .03, method="dual")  #1\% lower and 3\% upper trimmed least-
			           #squares fit. 
trq.print(trq(x,y)) #prints trq results in the style of ls.print.
}
\keyword{regression}
\section{METHOD}{
details of the methods may be found in Koenker and Portnoy(1987) for
the case of primal trimming and in Gutenbrunner and Jureckova(1991) for
dual trimming.
On the estimation of the covariance matrix for individual quantiles,
see Koenker(1987) and the discussion in Hendricks and Koenker(1991).
The estimation of the covariance matrix under  non-iid
conditions is an open research problem.  
}
\references{
Bassett, G., and Koenker, R. (1982), "An Empirical Quantile Function 
for Linear Models With iid Errors," 
\emph{Journal of the American Statistical Association,
}
77, 407-415.


Koenker, R.W. (1987), "A Comparison of Asymptotic Methods of Testing
based on L1 Estimation," in Y. Dodge (ed.) 
\emph{Statistical Data Analysis Based on the L1 norm and Related Methods,
}
New York:  North-Holland.


Koenker, R. W., and Bassett, G.W (1978), "Regression Quantiles",
\emph{Econometrica,
}
46, 33-50.


Koenker, R., and Portnoy, S. (1987), "L-Estimation for Linear Models",
\emph{Journal of the American Statistical Association,
}
82, 851-857.


Ruppert, D. and Carroll, R.J. (1980), "Trimmed Least Squares Estimation
in the Linear Model", 
\emph{Journal of the American Statistical Association, 75, 828-838.
}
}
\section{SEE ALSO}{
rq and qrq for further details.
}
% Converted by Sd2Rd version 0.3-2.
