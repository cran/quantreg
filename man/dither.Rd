
\name{dither}
\alias{dither}
\title{ Function to randomly perturb a vector}
\description{
With malice aforethought, dither adds a specified random perturbation to each element
of the input vector, usually employed as a device to mitigate the effect of ties.
}

\usage{
dither(x, type = "symmetric", value = NULL)
}
\arguments{
  \item{x}{\code{x} a numeric vector }
  \item{type}{\code{type} is either 'symmetric' or 'right' }
  \item{value}{\code{value} scale of dequantization }
}
\details{
The function \code{dither} operates slightly differently than the function
\code{jitter} in base R, permitting strictly positive perturbations with
the option \code{type = "right"} and using somewhat different default schemes
for the scale of the perturbation.  Dithering the response variable is
frequently a useful option in quantile regression fitting to avoid deleterious
effects of degenerate solutions.  See, e.g. Machado  and Santos Silva (2005).
For a general introduction and some etymology see the Wikipedia article on "dither".
For integer data it is usually advisable to use \code{value = 1}.
When 'x' is a matrix or array dither treats all elements as a vector but returns
an object of the original class.
}
\value{
A dithered version of the input vector 'x'.
}

\references{ 
Machado, J.A.F. and Santos Silva, J.M.C. (2005), Quantiles for Counts,  Journal of the American Statistical Association, vol. 100, no. 472, pp. 1226-1237. 
}
\author{ R. Koenker }
\note{ Some further generality might be nice, for example something other than
uniform noise would be desirable in some circumstances.  Note that when dithering
you are entering into the "state of sin" that John von Neumann famously attributed
to anyone considering "arithmetical methods of producing random digits."  If you 
need to preserve reproducibility, then \code{set.seed} is your friend.
}
\seealso{ \code{\link{jitter}} }
\examples{
x <- rlnorm(40)
y <- rpois(40, exp(.5 + log(x)))
f <- rq(dither(y, type = "right", value = 1) ~ x)
}
\keyword{ manip }

