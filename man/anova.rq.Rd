\name{anova.rq}
\alias{anova.rq}
\alias{anova.rqlist}
\alias{rq.test.rank}
\title{ Anova function for quantile regression fits }
\description{
	Compute test statistics for two or more quantile regression fits.
}
\usage{
anova.rq(object, ...)
anova.rqlist(object, ...,  test = "Wald", score = "tau")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object, ...}{ objects of class `rq', originating from a call to `rq'. }
  \item{test}{ A character string specifying the test statistic to use.
	Can be either `Wald' or `rank'.}
  \item{score}{ A character string specifying the score function to use,
	only needed or applicable for the `rank' form of the test.}
}
\details{
There are two (as yet) distinct forms of the test.  In the first the
fitted objects all have the same specified quantile (tau) and the intent
is to test the hypothesis that smaller models are adaquete relative to
the largest specified model.  In this case there are two options for
the argument `test', by default a Wald test is computed as in 
Bassett and Koenker (1982).  If test=`rank' is specified, then a rank
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
options.

The Wald form of the test is local in sense that the null hypothesis
asserts only that a subset of the covariates are ``insignificant'' at
the specified quantile of interest.  The rank form of the test can also be
used to test the global hypothesis that a subset is ``insignificant''
over an entire range of quantiles.  The use of the score function
score = "tau" restricts the rank test to the local hypothesis of
the Wald test.  

In the second form of the test the linear predictor of the fits are
all the same, but the specified quantiles (taus) are different.  In
this case the hypothesis of interest is that the slope coefficients of
the models are identical.  The test statistic is a variant of the Wald 
test described in Koenker and Bassett (1982).

By default, the tests return an F-like statistic in the sense that the
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
several thousand observations.
} 

\seealso{ The model fitting function  \code{\link{rq}}, 
and the functions for testing hypothesis on the entire quantile
regression process \code{\link{khmaladze.test}}.  For further details
on the rank tests see \code{\link{ranks}}}. 

\examples{
data(barro)
fit0 <- rq(y.net ~  lgdp2 + fse2 + gedy2 , data = barro)
fit1 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro)
fit2 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,tau=.75)
fit3 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,tau=.25)
anova(fit1,fit0)
anova(fit1,fit2,fit3)
}
\keyword{ htest }
\keyword{ regression}
\keyword{ robust}
