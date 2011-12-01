\name{lprq}
\alias{lprq}
\title{ locally polynomial quantile regression }
\description{
This is a toy function to illustrate how to do locally polynomial
quantile regression univariate smoothing.  
}
\usage{
lprq(x, y, h, tau = .5, m = 50)
}
\arguments{
  \item{x}{ The conditioning covariate}
  \item{y}{ The response variable  }
  \item{h}{ The bandwidth parameter }
  \item{tau}{ The quantile to be estimated }
  \item{m}{ The number of points at which the function is to be estimated }
}
\details{
The function obviously only does locally linear fitting but can be easily
adapted to locally polynomial fitting of higher order.  The author doesn't
really approve of this sort of smoothing, being more of a spline person,
so the code is left is its (almost) most trivial form. 
}
\value{
The function compute a locally weighted linear quantile regression fit
at each of the m design points, and returns:
  \item{xx}{The design points at which the evaluation occurs}
  \item{fv}{The estimated function values at these design points}
  \item{dev}{The estimated first derivative values at the design points}
}
\references{ Koenker, R. (2004) Quantile Regression }
\author{R. Koenker }
\note{One can also consider using B-spline expansions see \code{bs}.}

\seealso{ \code{rqss} for a general approach to oonparametric QR fitting.  }
\examples{
require(MASS)
data(mcycle)
attach(mcycle)
plot(times,accel,xlab = "milliseconds", ylab = "acceleration (in g)")
hs <- c(1,2,3,4)
for(i in hs){
        h = hs[i]
        fit <- lprq(times,accel,h=h,tau=.5)
        lines(fit$xx,fit$fv,lty=i)
        }
legend(50,-70,c("h=1","h=2","h=3","h=4"),lty=1:length(hs))
}
\keyword{smooth}
\keyword{robust}
