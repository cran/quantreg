\name{qss}
\alias{qss}
\alias{qss1}
\alias{qts1}
\alias{qss2}
\alias{triogram.fidelity}
\alias{triogram.penalty}
\title{Additive Nonparametric Terms for rqss Fitting}
\description{
In the formula specification of \code{rqss} nonparametric terms
are specified with \code{qss}.  Both univariate and bivariate
specifications are possible, and qualitative constraints may also be specified
for the qss terms.
}
\usage{
qss(x, constraint = "N", lambda = 1, ndum = 0, dummies = NULL, 
    Dorder = 1, w = rep(1, length(x)))
}
\arguments{
  \item{x}{The covariate determining the nonparametric component, if
	x is a matrix with two columns then  the qss
	function will construct a penalized triogram term.}
  \item{lambda}{The smoothing parameter governing the tradeoff between
	fidelity and the penalty component for this term.  Larger lambdas
	produce smoother fits. In future
	versions there should be an automatic mechanism for default
	choice of the lambdas.  For now, this is the responsibility
	of the user.}
  \item{constraint}{Optional specification of qualitative constraints
	on the fitted univariate qss functions, take the values: "N","I","D","V","C"
	"VI","VD","CI","CD" for none, increasing, decreasing, convex,
	concave, convex and increasing, etc.  And for bivariate qss
	components can take the values "N","V","C" for none, convex, and concave.
	Note that confidence bands for constrained fits of this sort, while
	available from \code{plot.rqss} as of yet lack a formal justification.}
  \item{ndum}{number of dummy vertices: this is only relevant for qss2
	terms.  In addition to vertices at the observed (x,y) points
	ndum dummy vertices are generated -- distributed uniformly over
	the rectangle given by the Cartesian product of the ranges of
	x and y -- observations that fall in the convex hull of the
	observations are retained. So the actual number of dummy
	vertices used is smaller than ndum.  The values of these
	vertices are returned in the list dummies, so that they can
	be reused.}
  \item{Dorder}{Order of the total variation penalty, the default of 1
      implies a penalty on the first derivative of the fitted function,
      a value of 0 implies total variation of the fitted function
      itself will be penalized.  Note that only monotonicity constraints,
      "I" and "D" are allowed when \code{Dorder = 0}, and result in estimates
      that are equivalent to a form of isotonic regression when lambda is
      sufficiently near zero.  Results in this case from the package \pkg{isotone}
      may differ slightly when plotted due to multiple solutions so it is prudent
      to evaluate the objective function for both solutions.}
  \item{dummies}{list of dummy vertices as generated, for example by
	triogram.fidelity when ndum > 0.  Should be a list with x
	and y components.  These points should lie inside the convex
	hull of the real xy points, but no explicit checking of this
	assertion is currently done.}
  \item{w}{ weights not yet unimplemented }
}
\details{
The various pieces returned are stored in sparse matrix.csr form.
See \code{\link{rqss}} for details on how they are assembled.  To preserve the
sparsity of the design matrix the first column of each qss term is dropped.
This differs from the usual convention that would have forced qss terms
to have mean zero.  This convention has implications for prediction that need to be
recognized.  The penalty components for qss terms are based on total
variation penalization of the first derivative (and gradient, for bivariate x)
as described in the references appearing in the help for \code{rqss}.
When Dorder = 0, fitting is like the taut string methods of Davies (2014), except
for the fact that fidelity is quantilesque rather than quadratic,
and that no provision is made for automatic selection of the smoothing 
parameter.

For the bivariate case, package \pkg{interp} (and for plotting also
\pkg{interp}) are required (automatically, by the \R code).
}
\value{
  \item{F}{Fidelity component of the design matrix}
  \item{dummies}{List of dummy vertices}
  \item{A}{Penalty component of the design matrix}
  \item{R}{Constraint component of the design matrix}
  \item{r}{Constraint component of the rhs}
}
\author{Roger Koenker}

\references{
    Davies, Laurie (2014) \emph{Data Analysis and Approximate Models}, CRC Press.
}
\seealso{ \code{\link{rqss}}}
\keyword{smooth}
\keyword{robust}
