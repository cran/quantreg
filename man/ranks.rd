\name{ranks}
\alias{ranks}
\title{
Function to compute ranks from the dual (regression rankscore) process
}
\description{
}
\usage{
ranks(v, score="wilcoxon")
}
\arguments{
\item{v}{
regression quantile structure for the model of interest
}
\item{score}{
The score function desired.  Currently implemented score functions are
Wilcoxon, Normal, and Sign which are asymptotically optimal for the
logistic, Gaussian and Laplace error models respectively.
}}
\value{
The function returns two components one is the ranks, the other is a
scale factor which is the L_2 norm of the score function.  All score
functions should be normalized to have mean zero.
}
\section{Side Effects}{
}
\details{
}
\references{
Gutenbrunner, C., J. Jureckova,  Koenker, R. and  Portnoy, S.(1993)
Tests of Linear Hypotheses  based on Regression Rank Scores",
Journal of Nonparametric Statistics, (2), 307-331.
}
\seealso{
See also rq, rrs.test
}
\examples{
data(stackloss)
ranks(rq(stack.x,stack.loss))
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-2.
