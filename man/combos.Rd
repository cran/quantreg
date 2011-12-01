\name{combos}
\alias{combos}
\title{Ordered Combinations}
\description{
All m combinations of the first n integers taken p at a time
are computed and return as an p by m matrix.  The columns
of the matrix are ordered so that adjacent columns differ
by only one element.  This is just a reordered version of
\code{combn} in base R, but the ordering is useful for some
applications.
}
\usage{
combos(n,p)
}
\arguments{
  \item{n}{The n in n choose p}
  \item{p}{The p in n choose p}
}
\value{
  a \code{matrix} of dimension p by \code{choose(n,p)}
} 
\note{
Implementation based on a Pascal algorithm of Limin Xiang 
and Kazuo Ushijima (2001) translated to ratfor for R.
If you have \pkg{rgl} installed you might try \code{demo("combos")}
for a visual impression of how this works.
}
\references{
Limin Xiang and Kazuo Ushijima (2001) 
"On O(1) Time Algorithms for Combinatorial Generation," 
\emph{Computer Journal}, 44(4), 292-302. 
}
\examples{
H <- combos(20,3)
}
\keyword{utilities}
