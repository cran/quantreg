\name{rqProcess}
\alias{rqProcess}
\title{ Compute Quantile Regression Process }
\description{
Computes the quantile regression process for the model
specified by the formula on the grid specified by the
taus argument.  Intended for use in \code{\link{khmaladze.test}}.
}
\usage{
rqProcess(formula, data, taus = seq(0.2, 0.8, by = 0.1))
}
\arguments{
  \item{formula}{ model formula  }
  \item{data}{data frame to be used to interpret formula }
  \item{taus}{ quantiles at which the process is to be evaluated}
}
\details{
The process computes point estimates and also returns an array
with the J and Hinv matrices needed to compute standardizations.
}
\value{
  \item{coefs}{Point estimates of the coefs of the model}
  \item{J}{J matrix of the usual qr sandwich formula, without
	any tau(1-tau) factor }
  \item{Hinv}{Hinv matrices of the usual qr sandwich formula,
	this is an array of dimension (p,p,m) }
}
\author{ R. Koenker }


\seealso{\code{\link{table.rq}}, \code{\link{khmaladze.test}}}
\examples{ }
\keyword{ regression }
\keyword{ htest }
