\name{QTECox}
\alias{QTECox}
\title{Function to obtain QTE from a Cox model}
\description{Computes quantile treatment effects comparable to those of 
crq model from a coxph object.}
\usage{
QTECox(x, smooth = TRUE)
}
\arguments{
  \item{x}{An object of class coxph produced by \code{coxph}.}
  \item{smooth}{Logical indicator if TRUE (default) 
	then Cox survival function is smoothed.}
}
\details{ Estimates of the Cox QTE, \eqn{\frac{dQ(t|x)}{dx_{j}}}{(d/dx_j) Q( t | x ) } 
          at \eqn{x=\bar{x}}{x=xbar}, can be expressed as a function of t as follows:

\deqn{\frac{dQ(t|x)}{dx_{j}}=\frac{dt}{dx_{j}}\frac{dQ(t|x)}{dt}}{
(d/dx_j) Q( t | x ) = (d/dx_j)t * (d/dt) Q(t | x)}

The Cox survival function, \eqn{S(y|x)=\exp \{-H_{0}(y)\exp (b^{\prime
}x)\}}{S( y | x ) = exp{ - H_o(y) exp(b'x) }}

\deqn{\frac{dS(y|x)}{dx_{j}}=S(y|x)log \{S(y|x)\}b_{j}}{(d/dx_j) 
S( y | x ) = S( y | x ) log(S( y | x )) b_j}


where \eqn{\frac{dQ(t|x)}{dx_{j}}}{ (d/dt) Q(t | x) }
can be estimated by \eqn{\frac{\Delta (t)}{\Delta (S)}
(1-t)}{- (diff(t)/diff(S) (1-t)}
where $S$  and $t$ denote the \code{surv} and \code{time} components
of the  \code{survfit} object.
Note that since \eqn{t=1-S(y|x)}{t = 1 - S( y | x )}, the above is the 
value corresponding to the argument $(1-t)$; and furthermore

\deqn{\frac{dt}{dx_{j}}=-\frac{dS(y|x)}{dx_{j}}=-(1-t) log (1-t)b_{j}}{
   (d/dx_j)t = - (d/dx_j) S( y | x ) =  - (1-t) log(1-t) b_j}

Thus the QTE at the mean of x's is:

\deqn{(1-S)= \frac{\Delta (t)}{\Delta (S)}S ~log
(S)b_{j}}{(1 - S) = (diff(t)/diff(S) S log(S) b_j}


Since \eqn{\Delta S}{diff(S)} is negative and $log (S)$ is also negative
this has the same sign as \eqn{b_{j}}
The crq  model fits the usual AFT form  Surv(log(Time),Status), then 

\deqn{\frac{d log (Q(t|x))}{dx_{j}}=\frac{dQ(t|x)}{dx_{j}}/
Q(t|x)}{(d/dx_j) log(Q( t | x )) = (d/dx_j) Q( t | x ) / Q( t | x )}

This is the matrix form returned.
}

\value{
  \item{taus }{points of evaluation of the QTE.}
  \item{QTE}{matrix of QTEs, the ith column contains the QTE for the
	ith covariate effect.  Note that there is no intercept effect.
	see \code{plot.summary.crqs} for usage.}
}

\references{Koenker, R. and Geling, O. (2001). Reappraising Medfly 
longevity: a quantile regression survival analysis, J. Amer. Statist. 
Assoc., 96, 458-468} 

\author{Roger Koenker Stephen Portnoy & Tereza Neocleous}
\seealso{\code{\link{crq}}}

\keyword{survival}
