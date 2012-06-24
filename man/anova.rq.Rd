\name{anova.rq}
\alias{anova.rq}
\alias{anova.rqlist}
\alias{print.anova.rq}
\alias{rq.test.rank}
\alias{rq.test.anowar}
\title{ Anova function for quantile regression fits }
\description{
	Compute test statistics for two or more quantile regression fits.
}
\usage{
\method{anova}{rq}(object, ...)
\method{anova}{rqlist}(object, ...,  test = "Wald", joint = TRUE, 
	score = "tau", se = "nid", R = 200, trim = NULL)
rq.test.rank(x0, x1, y, v = NULL, score = "wilcoxon", weights = NULL, tau=.5,
        iid = TRUE, delta0 = rep(0,NCOL(x1)), omega = 1, trim = NULL, pvalue = "F")
rq.test.anowar(x0,x1,y,tau,R)
\method{print}{anova.rq}(x, ...)
}
\arguments{
  \item{object, ...}{ objects of class `rq', originating from a call to `rq'. }
  \item{test}{ A character string specifying the test statistic to use.
	Can be either `Wald' or `rank'.}
  \item{joint}{A logical flag indicating whether tests of equality of slopes
	should be done as joint tests on all slope parameters, or whether
	(when joint = FALSE) separate tests on each of the slope parameters
	should be reported.  This option applies only to the tests of 
	equality of slopes in the case that estimated models correspond
	to distinct taus.}
  \item{score}{ A character string specifying the score function to use,
	only needed or applicable for the `rank' form of the test.}

  \item{trim}{optional trimming proportion parameter(s)  -- only applicable for the
	Wilcoxon score function --  when one value is provided there is symmetric
	trimming of the score integral to the interval \code{(trim, 1-trim)}, when
	there are two values provided, then the trimming restricts the integration
	to \code{(trim[1], trim[2])}.}
  \item{x}{ objects of class `summary.rq', originating from a call to `summary'. }
  \item{x0}{design matrix for the null component of the rank and anowar tests. }
  \item{x1}{design matrix for the alternative component of the rank and anowar tests. }
  \item{y}{response vector for the alternative component of the rank and anowar tests. }
  \item{v}{optional rq process fit}
  \item{se}{method for computing standard errors, either "nid" or "ker"}
  \item{tau}{quantile of interest for quantile specific forms of testing.}
  \item{iid}{logical flag for quantile specific forms of testing, if TRUE the
	test presumes that the conditional densities take identical values,
	if it is FALSE then local densities are estimated and used, see Koenker(2005)
	p. 90.}
  \item{delta0}{vector of hypothetical parameter values under test, typically zeros
	but can be specified to be nonzero in cases where simulations are being used
	to evaluate the validity of the non-central chisquare theory of the test.} 
  \item{omega}{value to be used for the score and F dependent constant appearing
	in the non-centrality parameter,  this is only needed/useful when delta0
	is specified to be non-zero.  In the usual Wilcoxon (untrimmed) case this
	value is the integral the squared density.}
  \item{pvalue}{type of p-value to be used, by default a pseudo F-statistic is
	produced and the corresponding F p-value is computed, otherwise the
	more conventional chisquared p-values are reported.}
  \item{weights}{optional weight vector to be used for fitting.}
  \item{R}{ The number of resampling replications for the anowar form of the test, 
	used to estimate the reference distribution for the test statistic.}
}
\details{
There are two (as yet) distinct forms of the test.  In the first the
fitted objects all have the same specified quantile (tau) and the intent
is to test the hypothesis that smaller models are adequate relative to
the largest specified model.  In the second form of the test the linear predictor 
of the fits are all the same, but the specified quantiles (taus) are different.  

In the former case there are three options for
the argument `test', by default a Wald test is computed as in 
Bassett and Koenker (1982).  If \code{test = 'anowar'} is specified
then the test is based on the procedure suggested in Chen, Ying, Zhang
and Zhao (2008); the test is based on the difference in the QR objective
functions at the restricted and unrestricted models with a reference
distribution computed by simulation.  The p-value of this form of the
test is produced by fitting a density to the simulation values forming
the reference distribution using the \code{logspline} function from
the \pkg{logspline} package.  The acronym anowar stands for  analysis
of weighted absolute residuals.  If \code{test='rank'} is specified, then a rank
test statistic is computed as described in Gutenbrunner, Jureckova,
Koenker and Portnoy (1993).  In the latter case one can also specify
a form for the score function of the rank test, by default the Wilcoxon
score is used, the other options are score=`sign' for median (sign) scores,
or score=`normal' for normal (van der Waerden) scores.  A fourth option
is score=`tau' which is a generalization of median scores to an arbitrary
quantile, in this case the quantile is assumed to be the one associated
with the fitting of the specified objects.  The computing of
the rank form of the test is carried out in the \code{\link{rq.test.rank}}
function, see \code{\link{ranks}} for further details on the score function
options.  The Wald form of the test is local in sense that the null hypothesis
asserts only that a subset of the covariates are ``insignificant'' at
the specified quantile of interest.  The rank form of the test can also be
used to test the global hypothesis that a subset is ``insignificant''
over an entire range of quantiles.  The use of the score function
score = "tau" restricts the rank test to the local hypothesis of
the Wald test.  

In the latter case the hypothesis of interest is that the slope coefficients of
the models are identical.  The test statistic is a variant of the Wald 
test described in Koenker and Bassett (1982).

By default, both forms of the tests return an F-like statistic in the sense that the
an asymptotically Chi-squared statistic is divided by its degrees of
freedom and the reported p-value is computed for an F statistic based on 
the numerator degrees of freedom equal to the rank of the null hypothesis and
the denominator degrees of freedom is taken to be the sample size
minus the number of parameters of the maintained model.
}
\value{
 An object of class `"anova"' inheriting from class `"data.frame"'.
}
\references{
[1] Bassett, G. and R. Koenker  (1982). Tests of Linear Hypotheses
and L1 Estimation,
\emph{Econometrica}, \bold{50}, 1577--83.

[2] Koenker, R. W. and Bassett, G. W. (1982). Robust Tests for
Heteroscedasticity based on Regression Quantiles,
\emph{Econometrica}, \bold{50}, 43--61.
 
[3] Gutenbrunner, C., Jureckova, J., Koenker, R, and S. Portnoy  (1993).
Tests of Linear Hypotheses based on Regression Rank Scores,
\emph{J. of Nonparametric Statistics}, \bold{2}, 307--331.

[4] Chen, K. Z. Ying, H. Zhang, and L Zhao, (2008) Analysis of least absolute
deviations, Biometrika, 95, 107-122.

[5] Koenker, R. W. (2005).  Quantile Regression, Cambridge U. Press. 
}
\author{ Roger Koenker }

\section{WARNING }{
An attempt to verify that the models are nested in the first form
of the test is made, but this relies on checking set inclusion of
the list of variable names and is subject to obvious ambiguities
when variable names are generic.
The comparison between two or more models will only be valid if
they are fitted to the same dataset. This may be a problem if
there are missing values and R's default of `na.action = na.omit' is used. 
The rank version of the nested model tests involves computing the entire
regression quantile process using parametric linear programming and thus
can be rather slow and memory intensive on problems with more than 
several thousand observations.} 
\seealso{The model fitting function  \code{\link{rq}}, 
and the functions for testing hypothesis on the entire quantile
regression process \code{\link{KhmaladzeTest}}.  For further details
on the rank tests see \code{\link{ranks}}.}
\examples{
data(barro)
fit0 <- rq(y.net ~  lgdp2 + fse2 + gedy2 , data = barro)
fit1 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro)
fit2 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,tau=.75)
fit3 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,tau=.25)
anova(fit1,fit0)
anova(fit1,fit2,fit3)
anova(fit1,fit2,fit3,joint=FALSE)
}
\keyword{ htest }
\keyword{ regression}
\keyword{ robust}
