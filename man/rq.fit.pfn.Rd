\name{rq.fit.pfn}
\alias{rq.fit.pfn}
\title{ Preprocessing Algorithm for Quantile Regression}
\description{
A preprocessing algorithm for the Frisch Newton algorithm 
for quantile regression. This is one possible method for rq().}
\usage{
rq.fit.pfn(x, y, tau=0.5, Mm.factor=0.8, max.bad.fixup=3, eps=1e-06)
}
\arguments{
  \item{x}{design matrix usually supplied via rq() }
  \item{y}{ response vector usually supplied via rq() }
  \item{tau}{ quantile of interest }
  \item{Mm.factor}{  constant to determine sub sample size m}
  \item{max.bad.fixup}{ number of allowed mispredicted signs of residuals }
  \item{eps}{ convergence tolerance }
}
\details{
Preprocessing algorithm to reduce the effective sample size for QR
problems with (plausibly) iid samples.  The preprocessing relies
on subsampling of the original data, so situations in which the
observations are not plausibly iid, are likely to cause problems.
The tolerance eps may be relaxed somewhat.
}
\value{
Returns an object of type rq
}
\references{ Portnoy and Koenker, Statistical Science, (1997) 279-300}
\author{ Roger Koenker <rkoenker@uiuc.edu>}
\seealso{ \code{\link{rq}}}

\keyword{ regression }
