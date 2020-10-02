\name{ParetoTest}
\alias{ParetoTest}
\alias{Hill}
\alias{Hill.fit}
\alias{print.Hill}
\alias{summary.Hill}
\alias{print.summary.Hill}
\alias{Pickands.fit}
\alias{Pickands}
\alias{print.Pickands}
\alias{summary.Pickands}
\alias{print.summary.Pickands}
\alias{Pickands.fit}
\title{Estimation and Inference on the Pareto Tail Exponent for Linear Models}
\description{Estimation and inference about the tail behavior of the response in
    linear models are based on the adaptation of the univariate Hill (1975)
    and Pickands (1975) estimators for quantile regression by Chernozhukov,
    Fernandez-Val and Kaji (2018).}

\usage{
ParetoTest(formula, tau = 0.1, data = NULL, flavor = "Hill", m = 2, cicov = .9, ...) 
}
\arguments{
\item{formula}{a formula specifying the model to fit by \code{\link{rq}}}
\item{tau}{A threshold on which to base the estimation}
\item{data}{a data frame within which to interpret the formula}
\item{flavor}{Currently limited to either "Hill" or "Pickands"}
\item{m}{a tuning parameter for the Pickands method .} 
\item{cicov}{Desired coverage probability of confidence interval.} 
\item{...}{other arguments to be passed to \code{\link{summary.rq}}.
    by default the summary method is the usual xy bootstrap, with
    \code{B = 200} replications.}
}
\value{
an object of class ParetoTest  is returned containing:

\item{z}{ A named vector with components: the estimate, a bias
corrected estimate, a lower bound of the confidence interval,
an upper bound of the confidence interval, and a Bootstrap
Standard Error estimate.}

\item{tau}{
    The tau threshold used to compute the estimate
}
}

\examples{
n = 500
x = rnorm(n)
y = x + rt(n,2)
Z = ParetoTest(y ~ x, .9, flavor = "Pickands")
}

\keyword{htest}
\references{
    Chernozhukov, Victor, Ivan Fernandez-Val, and Tetsuya Kaji, (2018)
    Extremal Quantile Regression, in Handbook of Quantile Regression, 
    Eds. Roger Koenker, Victor Chernozhukov, Xuming He, Limin Peng,
    CRC Press.

    Hill, B. M. (1975). A simple general approach to inference about the tail of a distribution. 
    The Annals of Statistics 3(5), 1163-1174.

    Pickands, J. (1975). Statistical inference using extreme order statistics. 
    The Annals of Statistics 3(1), 119-131.
}
    
