\name{bandwidth.rq}
\alias{bandwidth.rq}
\title{ bandwidth selection for rq functions }
\description{
 function to compute bandwidth for sparsity estimation 
}
\usage{
bandwidth.rq(p, n, hs=TRUE, alpha=0.05)
}
\arguments{
  \item{p}{ quantile(s) of interest }
  \item{n}{ sample size }
  \item{hs}{ flag for hall-sheather method }
  \item{alpha}{ alpha level for intended confidence intervals }
}
\details{ If hs=TRUE (default) then the Hall-Sheather(1988) rule \eqn{O(n^{-1/3})}
is used, if hs=FALSE then the Bofinger \eqn{O(n^{-1/5})} is used.
}
\value{
returns a vector of bandwidths corresponding to the argument p.
}
\references{ Hall and Sheather(1988, JRSS(B)),Bofinger (1975, Aus. J. Stat)}
\author{ Roger Koenker rkoenker@uiuc.edu}
\keyword{ regression }
