\name{plot.table.rq}
\alias{plot.table.rq}
\title{
Plot Table of Quantile Regression Results
}
\description{
The function makes an array of p - plots based on an array produced
by \code{table.rq} of quantile regression results.  The plots each represent
one parameter of the model specified in the formula argument to
\code{table.rq}; the plots consist of the point estimates betahat(tau)
plotted against the taus specified in the \code{table.rq} command with
a confidence band as produced by \code{rq()}.
}
\usage{
plot.table.rq(x,nrow=3,ncol=2, alpha=.1, ...)
}
\arguments{
\item{x}{
  object of class \code{table.rq} containing the array to be plotted.
}
\item{nrow}{
  number of rows per page of plots
}
\item{ncol}{
  number of columns per page of plots
}
\item{alpha}{
  alpha level of the confidence band
}
\item{...}{
  optional arguments to plot
}
}
\section{Side Effects}{
  plots an array of figures on the current graphics device.
}
\details{
  See \code{rq()} and \code{rq.fit.br()} and \code{table.rq()} for
  further details on control of data.
  Obviously, further plotting parameters could be added in a more full-blown
  version.  This version is meant just to be illustrative.
}
\seealso{
  \code{\link{rq}}, \code{\link{rq.fit.br}}, \code{\link{table.rq}}
}
\keyword{hplot}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
