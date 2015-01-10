\name{qrisk}
\alias{qrisk}
\title{ Function to compute Choquet portfolio weights}
\description{
This function solves a weighted quantile regression problem to find the
optimal portfolio weights minimizing a Choquet risk criterion described
in Bassett, Koenker, and Kordas (2002).
}
\usage{
qrisk(x, alpha = c(0.1, 0.3), w = c(0.7, 0.3), mu = 0.07, 
      R = NULL, r = NULL, lambda = 10000)
}
\arguments{
  \item{x}{n by q matrix of historical or simulated asset returns }
  \item{alpha}{vector of alphas receiving positive weights in the Choquet criterion}
  \item{w}{weights associated with alpha in the Choquet criterion  }
  \item{mu}{targeted rate of return for the portfolio}
  \item{R}{matrix of constraints on the parameters of the quantile regression, see below}
  \item{r}{rhs vector of the constraints described by R}
  \item{lambda}{Lagrange multiplier associated with the constraints}
}
\details{
The function calls \code{rq.fit.hogg} which in turn calls the constrained Frisch
Newton algorithm.  The constraints Rb=r are intended to apply only to the slope
parameters, not the intercept parameters.  The user is completely responsible to
specify constraints that are consistent, ie that have at least one feasible point.
See examples for imposing non-negative portfolio weights.
}
\value{
  \item{pihat}{the optimal portfolio weights}
  \item{muhat }{the in-sample mean return of the optimal portfolio}
  \item{qrisk}{the in-sample Choquet risk of the optimal portfolio}
}
\references{
http://www.econ.uiuc.edu/~roger/research/risk/risk.html

Bassett, G., R. Koenker, G Kordas, (2004) Pessimistic Portfolio Allocation and 
Choquet Expected Utility, J. of Financial Econometrics,  2, 477-492.
}
\author{ R. Koenker }
\examples{
#Fig 1:  ... of Choquet paper
        mu1 <- .05; sig1 <- .02; mu2 <- .09; sig2 <- .07
        x <- -10:40/100
        u <- seq(min(c(x)),max(c(x)),length=100)
        f1 <- dnorm(u,mu1,sig1)
        F1 <- pnorm(u,mu1,sig1)
        f2 <- dchisq(3-sqrt(6)*(u-mu1)/sig1,3)*(sqrt(6)/sig1)
        F2 <- pchisq(3-sqrt(6)*(u-mu1)/sig1,3)
        f3 <- dnorm(u,mu2,sig2)
        F3 <- pnorm(u,mu2,sig2)
        f4 <- dchisq(3+sqrt(6)*(u-mu2)/sig2,3)*(sqrt(6)/sig2)
        F4 <- pchisq(3+sqrt(6)*(u-mu2)/sig2,3)
        plot(rep(u,4),c(f1,f2,f3,f4),type="n",xlab="return",ylab="density")
        lines(u,f1,lty=1,col="blue")
        lines(u,f2,lty=2,col="red")
        lines(u,f3,lty=3,col="green")
        lines(u,f4,lty=4,col="brown")
        legend(.25,25,paste("Asset ",1:4),lty=1:4,col=c("blue","red","green","brown"))
#Now generate random sample of returns from these four densities.
n <- 1000
if(TRUE){ #generate a new returns sample if TRUE
	x1 <- rnorm(n)
	x1 <- (x1-mean(x1))/sqrt(var(x1))
	x1 <- x1*sig1 + mu1
	x2 <- -rchisq(n,3)
	x2 <- (x2-mean(x2))/sqrt(var(x2))
	x2 <- x2*sig1 +mu1
	x3 <- rnorm(n)
	x3 <- (x3-mean(x3))/sqrt(var(x3))
	x3 <- x3*sig2 +mu2
	x4 <- rchisq(n,3)
	x4 <- (x4-mean(x4))/sqrt(var(x4))
	x4 <- x4*sig2 +mu2
	}
library(quantreg)
x <- cbind(x1,x2,x3,x4)
qfit <- qrisk(x)
sfit <- srisk(x)
# Try new distortion function
qfit1 <- qrisk(x,alpha = c(.05,.1), w = c(.9,.1),mu = 0.09)
# Constrain portfolio weights to be non-negative
qfit2 <- qrisk(x,alpha = c(.05,.1), w = c(.9,.1),mu = 0.09,
	       R = rbind(rep(-1,3), diag(3)), r = c(-1, rep(0,3)))
}
\keyword{regression}
\keyword{robust}
\seealso{\code{\link{rq.fit.hogg}}, \code{\link{srisk}}}
