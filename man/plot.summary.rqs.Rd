\name{plot.summary.rqs}
\alias{plot.summary.rqs}

\title{Visualizing sequences of quantile regression summaries}

\description{A sequence of coefficient estimates for quantile
  regressions with varying \code{tau} parameters is visualized
  along with associated confidence bands.}

\usage{
\method{plot}{summary.rqs}(x, parm = NULL, level = 0.9, ols = TRUE,
  mfrow = NULL, mar = NULL, ylim = NULL, main = NULL,
  col = gray(c(0, 0.75)), border = NULL, lcol = 2, lty = 1:2,
  cex = 0.5, pch = 20, type = "b", xlab = "", ylab = "", \dots)
}

\arguments{
  \item{x}{an object of class \code{"summary.rqs"} as produce by 
    applying the \code{summary} method to a \code{\link{rq}} object
    (with a vector of \code{tau} values).}
  \item{parm}{a specification of which parameters are to be plotted,
    either a vector of numbers or a vector of names.  By default, all
    parameters are considered.}
  \item{level}{Confidence level of bands.  When using
    the rank based confidence intervals in summary, which is the default
    method for sample sizes under 1000, you will need to control the level
    of the intervals by passing the parameter alpha to 
    \code{\link{summary.rq}}, prior to calling 
    \code{\link{plot.summary.rqs}}.  Note also that alpha = 1 - level.}
  \item{ols}{logical. Should a line for the OLS coefficient and their confidence
    bands (as estimated by \code{\link[stats]{lm}}) be added?}
  \item{mfrow, mar, ylim, main}{graphical parameters. Suitable defaults are chosen
    based on the coefficients to be visualized.}
  \item{col}{vector of color specification for \code{rq} coefficients
    and the associated confidence polygon.}
  \item{border}{color specification for the confidence polygon. By default,
    the second element of \code{col} is used.}
  \item{lcol, lty}{color and line type specification for OLS coefficients
    and their confidence bounds.}
  \item{cex, pch, type, xlab, ylab, \dots}{further graphical parameters
    passed to \code{\link[graphics]{points}}.}
}

\details{The \code{plot} method for \code{"summary.rqs"} objects visualizes
  the coefficients along with their confidence bands. The bands can be
  omitted by using the \code{plot} method for \code{"rqs"} objects directly.}

\seealso{\code{\link{rq}}, \code{\link{plot.rqs}}}

\value{An array with all coefficients visualized (and associated confidence bands)
  is returned invisibly.}

\examples{
## fit Engel models (in levels) for tau = 0.1, ..., 0.9
data("engel")
fm <- rq(foodexp ~ income, data = engel, tau = 1:9/10)
sfm <- summary(fm)

## visualizations
plot(sfm)
plot(sfm, parm = 2, mar = c(5.1, 4.1, 2.1, 2.1), main = "", xlab = "tau", 
  ylab = "income coefficient", cex = 1, pch = 19)
}

\keyword{hplot}
