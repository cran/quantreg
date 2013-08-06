\name{sfn.control}
\alias{sfn.control}
\title{Set Control Parameters for Sparse Fitting }
\description{
Auxiliary function for setting storage dimensions and other parameters rq.fit.sfn[c]
}
\usage{
sfn.control(nsubmax = NULL, tmpmax = NULL, nnzlmax = NULL, cachsz = 64, 
	small = 1e-06, maxiter = 100, warn.mesg = TRUE)
}
\arguments{
  \item{nsubmax}{
upper bound for dimension of lindx
}
  \item{tmpmax}{
upper bound for dimension of tmpvec
}
  \item{nnzlmax}{
upper bound for non-zero entries of L stored in lnz, including diagonal
}

  \item{cachsz}{
size of cache in kbytes on target machine
}
  \item{small}{
convergence tolerance for interior point algorithm
}
  \item{maxiter}{
maximal number of interior point iterations. 
}
  \item{warn.mesg}{
logical flag controlling printing of warnings. 
}
}
\details{
Sparse fitting requires a number of temporary storage arrays whose size depends
on problem specific features in somewhat mysterious ways, parameters controlling
these sizes and some other fitting aspects can be controlled by specifying elements
of this control object.
}
\value{
List with components named as the arguments given above.
}

\author{
Roger Koenker
}

\seealso{
See Also \code{\link{rq.fit.sfn}}
}
\keyword{ utilities }
