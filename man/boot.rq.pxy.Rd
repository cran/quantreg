\name{boot.rq.pxy}
\alias{boot.rq.pxy}
\title{
    Preprocessing bootstrap method
}
\description{
    Bootstrap method exploiting preprocessing strategy to reduce 
    computation time for large problem.
}
\usage{
boot.rq.pxy(x, y, s, tau = 0.5, coef, method = "fn", Mm.factor = 3)
}
\arguments{
  \item{x}{
      Design matrix
}
  \item{y}{
      response vector
}
  \item{s}{
      matrix of multinomial draws for xy bootstrap
}
  \item{tau}{
      quantile of interest
}
  \item{coef}{
      point estimate of fitted object
}
  \item{method}{
      fitting method for bootstrap
}
  \item{Mm.factor}{
      constant to determine initial sample size
}
}
\details{
    See references for further details.
}
\value{
    Returns matrix of bootstrap estimates.
}
\references{
    Chernozhukov, V.  I. Fernandez-Val and B. Melly,
	Fast Algorithms for the Quantile Regression Process, 2019,
	arXiv, 1909.05782,

     Portnoy, S.  and R. Koenker, The Gaussian Hare and the Laplacian
	Tortoise, Statistical Science, (1997) 279-300
}
\author{
    Blaise Melly and Roger Koenker
}
\seealso{
\code{\link{rq.fit.ppro}}
}
\keyword{bootstrap} 
