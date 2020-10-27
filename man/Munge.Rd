\name{Munge}
\alias{Munge}
\title{
    Munge rqss formula
}
\description{
    function to recursively substitute arguments into rqss formula 
}
\usage{
Munge(formula, ...)
}
\arguments{
  \item{formula}{
      A rqss formula
}
  \item{\dots}{
      Arguments to be substituted into formula
}
}
\details{
    Intended (originally) for use with \code{demo(MCV)}.
    Based on an R-help suggestion of Gabor Grothendieck.
}
\value{
    A new formula after substitution
}
\seealso{
\code{demo(MCV)}
}
\examples{
lams <- c(1.3, 3.3)
f <- y ~ qss(x, lambda = lams[1]) + qss(z, lambda = lams[2]) + s
ff <- Munge(f, lams = lams)
}
\keyword{~manip}
