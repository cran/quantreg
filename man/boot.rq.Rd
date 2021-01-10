\name{boot.rq}
\alias{boot.rq}
\alias{boot.rq.xy}
\alias{boot.rq.wxy}
\alias{boot.rq.pwy}
\alias{boot.rq.spwy}
\alias{boot.rq.mcmb}
\title{ Bootstrapping Quantile Regression}
\description{
These functions can be used to construct standard errors, confidence
intervals and tests of hypotheses regarding quantile regression models.
}
\usage{
boot.rq(x, y, tau = 0.5, R = 200, bsmethod = "xy", mofn = length(y), 
	coef = NULL, blbn = NULL, cluster = NULL, U = NULL,  ...)
}
\arguments{
  \item{x}{ The regression design matrix}
  \item{y}{ The regression response vector}
  \item{tau}{ The quantile of interest}
  \item{R}{ The number of bootstrap replications}
  \item{bsmethod}{ The method to be employed.  There are (as yet) five
	options:  method = "xy" uses the xy-pair method, and
	method = "pwy" uses the method of Parzen, Wei and Ying (1994)
	method = "mcmb" uses the Markov chain marginal bootstrap
	of He and Hu (2002) and Kocherginsky, He and Mu (2003).
	The fourth method = "wxy" uses the generalized bootstrap
	of Bose and Chatterjee (2003) with unit exponential weights,
	see also Chamberlain and Imbens (2003).  The fifth method
	"wild" uses the wild bootstrap method proposed by Feng, He and Hu (2011). }
  \item{mofn}{ optional argument for the bootstrap method "xy" that
	permits subsampling (m out of n) bootstrap.  Obviously mofn
	should be substantially larger than the column dimension of x,
	and should be less than the sample size.}
  \item{coef}{coefficients from initial fitted object}
  \item{blbn}{orginal sample size for the BLB model}
  \item{cluster}{If non-NULL this argument should specify cluster id
      numbers for each observation, in which case the clustered version of 
      the bootstrap based on the proposal of Hagemann (2017). If present
      \code{bsmethod} is set to set to "cluster". If this option is used
      and the fitting method for the original call was "sfn" then the
      bootstrapping will be carried out with the "sfn" as well.  This
      is usually substantially quicker than the older version which
      employed the "br" variant of the simplex method.  Use of "sfn"
      also applies to the "pwy" method when the original fitting 
      was done with "sfn".  Finally, if \code{na.action = "omit"} and
      \code{length(object$na.action) > 0} then these elements are also
      removed from the \code{cluster} variable.  Consequently, the 
      length of the \code{cluster} variable should always be the same
      as the length of the original response variable before any 
      \code{na.action} takes place.  }
  \item{U}{If non-NULL this argument should specify an array of indices
      or gradient evaluations to be used by the corresponding bootstrap
      method as specified by \code{bsmethod}.  This is NOT intended as
      a user specified input, instead it is specified in \code{summary.rqs}
      to ensure that bootstrap samples for multiple taus use the same
      realizations of the random sampling.}
  \item{...}{ Optional arguments to control bootstrapping}
}
\details{
Their are several refinements that are still unimplemented.  Percentile
methods should be incorporated, and extensions of the methods to be used 
in anova.rq should be made.  And more flexibility about what algorithm is
used would also be good. 
}
\value{
A list consisting of two elements:
  A matrix \code{B} of dimension R by p is returned with the R resampled
	estimates of the vector of quantile regression parameters. When
	mofn < n for the "xy" method this matrix has been deflated by
	the factor sqrt(m/n).
  A matrix \code{U} of sampled indices (for \code{bsmethod in c("xy", "wxy")}) 
	or gradient evaluations (for \code{bsmethod in c("pwy", "cluster")})
	used to generate the bootstrapped realization, and potentially reused
	for other \code{taus} when invoked from \code{summary.rqs}.
}
\references{ 
[1] Koenker, R. W. (1994). Confidence Intervals for regression quantiles, in
P. Mandl and M. Huskova (eds.), \emph{Asymptotic Statistics}, 349--359,
Springer-Verlag, New York.

[2] Kocherginsky, M., He, X. and Mu, Y. (2005).
Practical Confidence Intervals for Regression Quantiles,
Journal of Computational and Graphical Statistics, 14, 41-55.

[3] Hagemann, A. (2017) Cluster Robust Bootstrap inference in 
quantile regression models, Journal of the American Statistical Association , 
112, 446--456.

[4] He, X. and Hu, F. (2002). Markov Chain Marginal Bootstrap.
Journal of the American Statistical Association , Vol. 97, no. 459,
783-795. 

[5] Parzen, M. I., L. Wei,  and Z. Ying  (1994): A resampling
method based on pivotal estimating functions,'' Biometrika, 81, 341--350.

[6] Bose, A. and S. Chatterjee, (2003) Generalized bootstrap for estimators
of minimizers of convex functions, \emph{J. Stat. Planning and Inf}, 117, 225-239.

[7]  Chamberlain G.  and Imbens G.W.  (2003) Nonparametric Applications of 
Bayesian Inference, Journal of Business & Economic Statistics, 21, pp. 12-18.

[8]  Feng, Xingdong, Xuming He, and Jianhua Hu (2011) Wild Bootstrap for
Quantile Regression, Biometrika, 98, 995--999. 
}

\author{ Roger Koenker (and Xuming He and M. Kocherginsky for the mcmb code)}
\seealso{  \code{\link{summary.rq}}}
\examples{
y <- rnorm(50)
x <- matrix(rnorm(100),50)
fit <- rq(y~x,tau = .4)
summary(fit,se = "boot", bsmethod= "xy")
summary(fit,se = "boot", bsmethod= "pwy")
#summary(fit,se = "boot", bsmethod= "mcmb")
}
\keyword{ regression}
