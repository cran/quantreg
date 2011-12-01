
\name{critval}
\alias{critval}
\title{
Hotelling Critical Values
}
\description{
Critical values for uniform confidence bands for rqss fitting
}
\usage{
critval(kappa, alpha = 0.05, rdf = 0)
}
\arguments{
\item{kappa}{
  length of the tube
}
\item{alpha}{
  desired non-coverage of the band, intended coverage is 1 - alpha
}
\item{rdf}{
  "residual" degrees of freedom of the fitted object.  If \code{rdf=0}
  then the Gaussian version of the critical value is computed, otherwise
  the value is based on standard Student t theory.
}
}
\value{
  A scalar critical value that acts as a multiplier for the uniform
  confidence band construction.
}
\details{
  The Hotelling tube approach to inference has a long and illustrious
  history.  See Johansen and Johnstone (1989) for an overview.  The implementation
  here is based on Sun and Loader (1994) and Loader's \pkg{locfit} package, although
  a simpler root finding approach is substituted for the iterative method used
  there.  At this stage, only univariate bands may be constructed.
}
\references{
Hotelling, H.  (1939): ``Tubes and Spheres in $n$-spaces, and a class
  of statistical problems,'' \emph{Am J. Math}, 61, 440--460.

Johansen, S.,   I.M. Johnstone  (1990): ``Hotelling's
  Theorem on the Volume of Tubes: Some Illustrations in Simultaneous
  Inference and Data Analysis,'' \emph{The Annals of Statistics}, 18, 652--684.

Sun, J. and C.V. Loader:  (1994) ``Simultaneous Confidence Bands for Linear Regression
and smoothing,'' \emph{The Annals of Statistics}, 22, 1328--1345.
}
\seealso{
  \code{\link{plot.rqss}}
}
\keyword{regression}
