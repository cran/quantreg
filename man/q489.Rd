\name{q489}
\alias{q489}
\title{Even Quicker Sample Quantiles }
\description{
The function \code{q489} computes a single sample quantile using a 
fortran implementation of the Floyd and Rivest (1975) algorithm.
In contrast to the more elaborate function \code{kuantile} that uses
the Kiweil (2005) implementation it does not attempt to replicate the
nine varieties of quantiles as documented in the base function.
\code{quantile}
}
\usage{
q489(x, tau = .5) 
}
\arguments{
  \item{x}{numeric vector}
  \item{tau}{the quantile of intereste.}
}
\details{ This is a direct translation of the Algol 68 implementation of
Floyd and Rivest (1975), implemented in Ratfor.  For the median, average 
case behavior requires \eqn{1.5 n + O((n log n)^{1/2})} comparisons. 
In preliminary experiments it seems to be somewhat faster in large samples
than the implementation \code{kuantile} of Kiwiel (2005). See  Knuth (1998)
for further details. No provision is made for non-uniqueness of the quantile.
so, when \eqn{\tau n} is an integer there may be some discrepancy.}  
\value{
  A scalar quantile of the same length as the vector p.
}
\references{ 
R.W. Floyd and R.L. Rivest: "Algorithm 489: The Algorithm
        SELECT---for Finding the $i$th Smallest of $n$ Elements",
        Comm. ACM 18, 3 (1975) 173,

K.C. Kiwiel: On Floyd and Rivest's SELECT Algorithm, Theoretical
      Computer Sci. 347 (2005) 214-238.

D. Knuth, The Art of Computer Programming, Volume 3, Sorting and 
	Searching, 2nd Ed., (1998), Addison-Wesley.
}
\author{ R.W.Floyd and R.L.Rivest, R implementation:  Roger Koenker }
\seealso{\code{\link{quantile}}}
\examples{
     medx <- q489(rnorm(1001))
}
\keyword{univar}
