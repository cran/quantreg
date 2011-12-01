\name{boot.rq}
\alias{boot.rq}
\alias{boot.rq.xy}
\alias{boot.rq.wxy}
\alias{boot.rq.pwy}
\alias{boot.rq.mcmb}
\title{ Bootstrapping Quantile Regression}
\description{
These functions can be used to construct standard errors, confidence
intervals and tests of hypotheses regarding quantile regression models.
}
\usage{
boot.rq(x, y, tau = 0.5, R = 200, bsmethod = "xy", mofn = length(y), ...)
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
  \item{...}{ Optional arguments to control bootstrapping}
}
\details{
Their are several refinements that are still unimplemented.  Percentile
methods should be incorporated, and extensions of the methods to be used 
in anova.rq should be made.
}
\value{
  A matrix of dimension R by p is returned with the R resampled
  estimates of the vector of quantile regression parameters. When
  mofn < n for the "xy" method this matrix has been deflated by
  the factor sqrt(m/n)
}
\references{ 

[1] Koenker, R. W. (1994). Confidence Intervals for regression quantiles, in
P. Mandl and M. Huskova (eds.), \emph{Asymptotic Statistics}, 349--359,
Springer-Verlag, New York.

[2] Kocherginsky, M., He, X. and Mu, Y. (2005).
Practical Confidence Intervals for Regression Quantiles,
Journal of Computational and Graphical Statistics, 14, 41-55.

[3] He, X. and Hu, F. (2002). Markov Chain Marginal Bootstrap.
Journal of the American Statistical Association , Vol. 97, no. 459,
783-795. 

[4] Parzen, M. I., L. Wei,  and Z. Ying  (1994): A resampling
method based on pivotal estimating functions,'' Biometrika, 81, 341--350.

[5] Bose, A. and S. Chatterjee, (2003) Generalized bootstrap for estimators
of minimizers of convex functions, \emph{J. Stat. Planning and Inf}, 117, 225-239.

[6]  Chamberlain G.  and Imbens G.W.  (2003) Nonparametric Applications of 
Bayesian Inference, Journal of Business & Economic Statistics, 21, pp. 12-18.

[7]  Feng, Xingdong, Xuming He, and Jianhua Hu (2011) Wild Bootstrap for
Quantile Regression, Biometrika, to appear.
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
