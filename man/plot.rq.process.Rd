\name{plot.rq}
\alias{plot.summary.rqs}
\alias{plot.rq.process}
\title{ plot the coordinates of the quantile regression process}
\description{ Function to plot quantile regression  process.  }
\usage{plot.summary.rqs(x, nrow=3, ncol=2, alpha = 0.1, ols = TRUE, ...) }
\arguments{
  \item{x}{ an object produced by rq() fitting  }
  \item{nrow}{ rows in mfrow }
  \item{ncol}{ columns in mfrow}
  \item{alpha}{alpha level of confidence bands}
  \item{ols}{if TRUE, plot ols coefs and their error band}
  \item{...}{ optional arguments to plot}
}
\author{ Roger Koenker rkoenker@uiuc.edu}
\seealso{ \code{\link{rq}}}
\keyword{hplot}
