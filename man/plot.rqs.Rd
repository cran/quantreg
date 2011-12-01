\name{plot.rqs}
\alias{plot.rqs}

\title{Visualizing sequences of quantile regressions}

\description{A sequence of coefficient estimates for quantile
  regressions with varying \code{tau} parameters is visualized.}

\usage{
\method{plot}{rqs}(x, parm = NULL, ols = TRUE,
  mfrow = NULL, mar = NULL, ylim = NULL, main = NULL, col = 1:2, lty = 1:2,
  cex = 0.5, pch = 20, type = "b", xlab = "", ylab = "", \dots)
}

\arguments{
  \item{x}{an object of class \code{"rqs"} as produce by \code{\link{rq}}
    (with a vector of \code{tau} values).}
  \item{parm}{a specification of which parameters are to be plotted,
    either a vector of numbers or a vector of names.  By default, all
    parameters are considered.}
  \item{ols}{logical. Should a line for the OLS coefficient (as estimated
    by \code{\link[stats]{lm}}) be added?}
  \item{mfrow, mar, ylim, main}{graphical parameters. Suitable defaults are chosen
    based on the coefficients to be visualized.}
  \item{col, lty}{graphical parameters. For each parameter, the first
    element corresponds to the \code{rq} coefficients and the second to
    the \code{lm} coefficients.}
  \item{cex, pch, type, xlab, ylab, \dots}{further graphical parameters
    passed to \code{\link[graphics]{plot}}.}
}

\details{The \code{plot} method for \code{"rqs"} objects visualizes the
  coefficients only, confidence bands can be added by using the \code{plot}
  method for the associated \code{"summary.rqs"} object.}

\seealso{\code{\link{rq}}, \code{\link{plot.summary.rqs}}}

\value{A matrix with all coefficients visualized is returned invisibly.}

\examples{
## fit Engel models (in levels) for tau = 0.1, ..., 0.9
data("engel")
fm <- rq(foodexp ~ income, data = engel, tau = 1:9/10)

## visualizations
plot(fm)
plot(fm, parm = 2, mar = c(5.1, 4.1, 2.1, 2.1), main = "", xlab = "tau", 
  ylab = "income coefficient", cex = 1, pch = 19)
}

\keyword{hplot}
