\name{rq.fit.br}
\alias{rq.fit.br}
\title{
Quantile Regression Fitting by Exterior Point Methods
}
\description{
  This function controls the details of QR fitting by the simplex approach
  embodied in the algorithm of Koenker and d'Orey based on the median
  regression algorithm of Barrodale and Roberts.  Typically, options
  controlling the construction of the confidence intervals would be passed
  via the \code{\dots{}} argument of \code{rq()}.
}
\usage{
rq.fit.br(x, y, tau=0.5, alpha=0.1, ci=FALSE, iid=TRUE, interp=TRUE, tcrit=TRUE)
}
\arguments{
\item{x}{
  the design matrix
}
\item{y}{
  the response variable
}
\item{tau}{
  the quantile desired, if tau lies outside (0,1) the whole process
  is estimated.
}
\item{alpha}{
  the nominal noncoverage probability for the confidence intervals, i.e. 1-alpha
  is the nominal coverage probability of the intervals. 
}
\item{ci}{
  logical flag if T then compute confidence intervals for the parameters
  using the rank inversion method of Koenker (1994).  See \code{rq()} for more
  details.  If F then return only the estimated coefficients.  Note that
  for large problems the default option ci = TRUE can be rather slow.
  Note also that rank inversion only works for p>1, an error message is
  printed in the case that ci=T and p=1.
}
\item{iid}{
  logical flag if T then the rank inversion is based on an assumption of
  iid error model, if F then it is based on an nid error assumption.
  See Koenker and Machado (1999) for further details on this distinction.
}
\item{interp}{
  As with typical order statistic type confidence intervals the test
  statistic is discrete, so it is reasonable to consider intervals that
  interpolate between values of the parameter just below the specified
  cutoff and values just above the specified cutoff.  If \code{interp =
    F} then
  the 2 ``exact'' values above and below on which the interpolation would
  be based are returned.
}
\item{tcrit}{
Logical flag if T -  Student t critical values are used, if F then normal
values are used.
}
}
\value{
  Returns an object of class \code{"rq"}
  for tau in (0,1), or else of class \code{"rq.process"}.
  Note that \code{rq.fit.br} when called for a single tau value
  will return the vector of optimal dual variables.
  See \code{\link{rq.object}} and \code{\link{rq.process.object}}
  for further details.
}
\details{
  If tau lies in (0,1) then an object of class \code{"rq"} is
  returned with various
  related inference apparatus.  If tau lies outside [0,1] then an object
  of class \code{rq.process} is returned.  In this case parametric programming
  methods are used to find all of the solutions to the QR problem for
  tau in (0,1), the p-variate resulting process is then returned as the
  array sol containing the primal solution and dsol containing the dual
  solution.  There are roughly \eqn{O(n \log n))}{O(n log n)} distinct
  solutions, so users should
  be aware that these arrays may be large and somewhat time consuming to
  compute for large problems.
}
\references{
Koenker, R. and J.A.F. Machado, (1999) Goodness of fit and related inference 
processes for quantile regression,
\emph{J. of Am Stat. Assoc.}, 94, 1296-1310. 
}
\seealso{
\code{\link{rq}}, \code{\link{rq.fit.fnb}}
}
\examples{
data(stackloss)
rq.fit.br(stack.x, stack.loss, tau=.73 ,interp=FALSE)
}
\keyword{regression}
