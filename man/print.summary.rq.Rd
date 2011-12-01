\name{print.summary.rq}
\alias{print.summary.rq}
\alias{print.summary.rqs}
\title{ Print Quantile Regression Summary Object }
\usage{
\method{print}{summary.rq}(x, digits=max(5, .Options$digits - 2), ...)
\method{print}{summary.rqs}(x, ...)
}
\arguments{
  \item{x}{
    This is an object of class \code{"summary.rq"} produced by a call to 
    \code{summary.rq()}.
  }
  \item{digits}{
    Significant digits reported in the printed table.
  }
  \item{...}{
    Optional arguments passed to printing function
  }
}
\description{Print summary of quantile regression object}
\seealso{  \code{\link{summary.rq}} }

\keyword{ regression }
