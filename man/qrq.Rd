\name{qrq}
\alias{qrq}
\title{
Linearized Quantile Estimation
}
\description{
Compute linearized quantiles from rq data structure.
}
\usage{
qrq(s, a)
}
\arguments{
\item{s}{
data structure returned by the quantile regression function rq with t<0 or t>1.
}
\item{a}{
the vector of quantiles for which the corresponding linearized quantiles are to be computed.
}}
\value{
a vector of the linearized quantiles corresponding to vector a,
as interpolated from the second row of s$sol.
}
\seealso{
rq and  trq  for further detail.
}
\examples{
x <- -10:10; y <- 0.2 * x + rt(x, df=3)
z <- qrq(rq(x,y), a=0.1)       # assigns z the linearized quantiles 
                               # corresponding to vector a.
}
\keyword{}
% Converted by Sd2Rd version 0.3-2.
