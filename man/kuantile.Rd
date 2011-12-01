\name{kuantile}
\alias{kuantile}
\alias{kselect}
\alias{kunique}
\title{Quicker Sample Quantiles }
\description{
The function 'kuantile' computes sample quantiles corresponding
to the specified probabilities. The intent is to mimic the generic
(base) function 'quantile' but using a variant of the Floyd and
Rivest (1975) algorithm which is somewhat quicker, especially for
large sample sizes.
}
\usage{
kuantile(x, probs = seq(0, 1, .25), na.rm = FALSE, names = TRUE, type = 7, ...) 
}
\arguments{
  \item{x}{numeric vector whose sample quantiles are wanted.}
  \item{probs}{numeric vector of probabilities with values in [0,1].}
  \item{type}{ an integer between 1 and 9 selecting one of the nine quantile
          algorithms detailed below to be used.}
  \item{na.rm}{logical: if true, any 'NA' and 'NaN''s are removed from 'x'
          before the quantiles are computed.}
  \item{names}{logical: if true, the result has a 'names' attribute. } 
  \item{...}{further arguments passed to or from other methods.}
}
\details{ A vector of length 'length(p)' is returned.  See the documentation
for 'quantile' for further details on the types.  The algorithm was written
by K.C. Kiwiel.  It is a modified version of the (algol 68) SELECT procedure of
Floyd and Rivest (1975), incorporating modifications of Brown(1976).
The algorithm has linear growth in the number of comparisons required as
sample size grows.  For the median, average case behavior requires
\eqn{1.5 n + O((n log n)^{1/2})} comparisons. See Kiwiel (2005) and Knuth (1998)
for further details.  When the number of required elements of p is large, it
may be preferable to revert to a full  sort.}
\value{
  A vector of quantiles of the same length as the vector p.
}
\references{ 
R.W. Floyd and R.L. Rivest: "Algorithm 489: The Algorithm
        SELECT---for Finding the $i$th Smallest of $n$ Elements",
        Comm. ACM 18, 3 (1975) 173,

T. Brown: "Remark on Algorithm 489", ACM Trans. Math.
        Software 3, 2 (1976), 301-304.

K.C. Kiwiel: On Floyd and Rivest's SELECT Algorithm, Theoretical
      Computer Sci. 347 (2005) 214-238.

D. Knuth, The Art of Computer Programming, Volume 3, Sorting and 
	Searching, 2nd Ed., (1998), Addison-Wesley.
}
\author{ K.C. Kiwiel, R interface:  Roger Koenker }
\seealso{\code{\link{quantile}}}
\examples{
     kuantile(x <- rnorm(1001))# Extremes & Quartiles by default

     ### Compare different types
     p <- c(0.1,0.5,1,2,5,10,50)/100
     res <- matrix(as.numeric(NA), 9, 7)
     for(type in 1:9) res[type, ] <- y <- kuantile(x,  p, type=type)
     dimnames(res) <- list(1:9, names(y))
     ktiles <- res

     ### Compare different types
     p <- c(0.1,0.5,1,2,5,10,50)/100
     res <- matrix(as.numeric(NA), 9, 7)
     for(type in 1:9) res[type, ] <- y <- quantile(x,  p, type=type)
     dimnames(res) <- list(1:9, names(y))
     qtiles <- res

     max(abs(ktiles - qtiles))


}
\keyword{univar}
