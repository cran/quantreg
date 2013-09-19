\name{boot.crq}
\alias{boot.crq}
\title{ Bootstrapping Censored Quantile Regression}
\description{
Functions used to estimated standard errors, confidence
intervals and tests of hypotheses for censored quantile regression models
using the Portnoy and Peng-Huang methods.  }
\usage{
boot.crq(x, y, c, taus, method, ctype = "right", R = 100, mboot, bmethod = "jack", ...)
}
\arguments{
  \item{x}{ The regression design matrix}
  \item{y}{ The regression response vector}
  \item{c}{ The censoring indicator}
  \item{taus}{ The quantiles of interest}
  \item{method}{ The fitting method: either "P" for Portnoy or "PH" for Peng and Huang.}
  \item{ctype}{ Either "right" or "left"}
  \item{R}{ The number of bootstrap replications}
  \item{bmethod}{ The bootstrap method to be employed.  There are (as yet) three
	options:  method  = "jack" uses the delete-d jackknife method
	described by Portnoy (2013), method = "xy-pair" uses the xy-pair method, 
	that is the usual multinomial resampling of xy-pairs, while  method
	= "Bose" uses the Bose and Chatterjee (2003) weighted resampling
	method with exponential weights.  The "jack" method is now the default.}
  \item{mboot}{ optional argument for the bootstrap method:  for bmethod = "jack"
	 it specifies the number, d, of the delete-d jackknife, for 
	 method = "xy-pair" it specifies the size of the bootstrap samples,
	 that permits subsampling (m out of n) bootstrap.  By default in the
	 former case it is set to 2 [sqrt(n)], for the latter the default is
	 n.  Obviously mboot should be substantially larger than the column dimension 
	 of x, and should be less than the sample size in both cases.}
  \item{...}{ Optional further arguments to control bootstrapping}
}
\details{
There are several refinements that are still unimplemented.  Percentile
methods should be incorporated, and extensions of the methods to be used 
in anova.rq should be made.  Note that bootstrapping for the Powell 
method "FP" is done via \code{\link{boot.rq}}.
}
\value{
  A matrix of dimension R by p is returned with the R resampled
  estimates of the vector of quantile regression parameters. When
  mofn < n for the "xy" method this matrix has been deflated by
  the fact sqrt(m/n)
}
\references{ 
 Bose, A. and S. Chatterjee, (2003) Generalized bootstrap for estimators
        of minimizers of convex functions, \emph{J. Stat. Planning and Inf}, 117,
        225-239.
 Portnoy, S. (2013) The Jackknife's Edge:  Inference for Censored Quantile Regression,
	\emph{CSDA}, forthcoming.

}

\author{ Roger Koenker }
\seealso{  \code{\link{summary.crq}}}
\keyword{ regression}
