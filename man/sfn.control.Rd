\name{sfn.control}
\alias{sfn.control}
\title{Set Control Parameters for Sparse Fitting}
\description{
  Auxiliary function for setting storage dimensions and other parameters
  for \code{\link{rq.fit.sfn}()}, \code{\link{rq.fit.sfnc}()} and \code{\link{rqss}()}.
}
\usage{
sfn.control(nsubmax = NULL, tmpmax = NULL, nnzlmax = NULL, cachsz = 64,
            small = 1e-06, maxiter = 100,
            tiny = 1e-30, Large = 1e128,
            warn.mesg = TRUE)
}
\arguments{
  \item{nsubmax}{upper bound for dimension of lindx}
  \item{tmpmax}{upper bound for dimension of tmpvec}
  \item{nnzlmax}{upper bound for non-zero entries of L stored in lnz, including diagonal}
  \item{cachsz}{size of cache in kbytes on target machine}
  \item{small}{convergence tolerance for interior point algorithm}
  \item{maxiter}{maximal number of interior point iterations}
  \item{tiny}{a tiny positive number; values below \code{tiny * max(diag)} are replaced by \code{Large};
    originally was \eqn{10^{-30}} hardcoded in Fortran code.}
  \item{Large}{a large number, practically \dQuote{Infinite} to replace
    \code{tiny} diagonal entries in Cholesky; was \eqn{10^{128}}, hardcoded in
    compiled code.}
  \item{warn.mesg}{logical flag controlling printing of warnings.}
}
\details{
Sparse fitting requires a number of temporary storage arrays whose size depends
on problem specific features in somewhat mysterious ways, parameters controlling
these sizes and some other fitting aspects can be controlled by specifying elements
of this control object.
}
\value{
  A \code{\link{list}} with components named as the arguments given above.
}

\author{
Roger Koenker
}

\seealso{
  \code{\link{rq.fit.sfn}},
  \code{\link{rq.fit.sfnc}},
  and \code{\link{rqss}} from which \code{sfn.control()} is called.
}
\keyword{ utilities }
