\name{table.rq}
\alias{table.rq}
\title{
  Table of Quantile Regression Results
}
\description{
  Function to produce a table of quantile regression results for a group
  of specified quantiles.
}
\usage{
table.rq(formula, taus=c(0.25, 0.5, 0.75), method="br", \dots)
}
\arguments{
  \item{formula}{
    formula for the linear model, see \code{rq()}
  }
  \item{taus}{
    quantiles of interest
  }
  \item{method}{
    algorithmic method, for large problems \code{method="fn"}
    would be preferred.
  }
  \item{...}{
    other optional arguments passed to \code{rq()}.
  }
}
\value{
  The function returns an array with dimensions (p, m, 3), where p is the
  dimension of the parameter vector of the model specified by formula, m
  is the number of quantiles specified by tau.  For each coefficient at
  each tau there is a point estimate and lower and upper limits of a
  confidence interval.  The object returned is of class \code{table.rq}
  and can be plotted, or formatted into a latex style table.
}
\details{
  This is only implemented for \code{method="br"}, but modifications for
  \code{"fn"} would be quite straightforward.  There is also an
  implementation for making a latex table in the Splus version,
  but this wasn't incorporated into the \R package due to some
  incompatibilities in the unix/system commands.
}
\seealso{
\code{\link{rq}}, \code{\link{rq.fit.br}}, \code{\link{plot.table.rq}}
}
\examples{
data(stackloss)
plot(table.rq(stack.loss~stack.x))#plot results of a quantile regression
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.


