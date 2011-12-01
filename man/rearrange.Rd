\name{rearrange}
\alias{rearrange}
\title{Rearrangement}
\description{ Monotonize a step function by rearrangement }
\usage{ rearrange(f,xmin,xmax) }
\arguments{ 
   \item{f}{ object of class stepfun } 
   \item{xmin}{minimum of the support of the rearranged f} 
   \item{xmax}{maximum of the support of the rearranged f} 

}
\details{
   Given a stepfunction \eqn{Q(u)}, not necessarily monotone, let 
   \eqn{F(y) = \int \{ Q(u) \le y \} du} denote the associated cdf
   obtained by randomly evaluating \eqn{Q} at \eqn{U \sim U[0,1]}.  The
   rearranged version of \eqn{Q} is \eqn{\tilde Q (u) = \inf \{
   u: F(y) \ge u \}.  The rearranged function inherits the right
   or left continuity of original stepfunction.}
}
\value{ Produces transformed stepfunction that is monotonic increasing.  }
\references{ 
Chernozhukov, V., I. Fernandez-Val, and A. Galichon, (2006)  Quantile and Probability
Curves without Crossing, Econometrica, forthcoming.

Chernozhukov, V., I. Fernandez-Val, and A. Galichon, (2009)  Improving Estimates of
Monotone Functions by Rearrangement, Biometrika, 96, 559--575. 

Hardy, G.H., J.E. Littlewood, and G. Polya (1934)  Inequalities,  Cambridge U. Press.
}
\author{R. Koenker}
\seealso{ \code{\link{rq}} \code{\link{rearrange}}}
\examples{
data(engel)
z <- rq(foodexp ~ income, tau = -1,data =engel)
zp <- predict(z,newdata=list(income=quantile(engel$income,.03)),stepfun = TRUE)
plot(zp,do.points = FALSE, xlab = expression(tau),
        ylab = expression(Q ( tau )), main="Engel Food Expenditure Quantiles")
plot(rearrange(zp),do.points = FALSE, add=TRUE,col.h="red",col.v="red")
legend(.6,300,c("Before Rearrangement","After Rearrangement"),lty=1,col=c("black","red"))
}
\keyword{regression}
