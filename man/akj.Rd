\name{akj}
\alias{akj}
\title{
Density estimation using adaptive kernel method
}
\description{univariate adaptive kernel density estimation a la Silverman}
\usage{
akj(x, z, p, h, alpha, kappa, iker1, iker2)


}
\arguments{
\item{x}{
points used for centers of kernel assumed to be sorted}
\item{z}{
points at which density is calculated; default to seq( min(x), max(x), 2*length(x) )}
\item{p}{
vector of probabilities associated with x's; default to 1/len(x) for each x.}
\item{h}{
initial window size (overall); default to Silverman's normal reference}
\item{alpha}{a sensitivity parameter that determines the sensitivity of the local bandwidth to variations in the pilot density; default to .5}
\item{kappa}{
constant determining initial (default) window width}
\item{iker1}{
kernel indicator, 0 for normal kernel (default) while 1 for cauchy kernel}
\item{iker2}{ xxx}
}
\value{
a R structure is returned

\item{dens}{
the vector of estimated density}
\item{psi}{
a vector of \eqn{\psi=-f'/f} function}
\item{score}{
a vector of score \eqn{(f'/f)^2-f'/f} function}
\item{h}{
same as the input argument h}
}
\keyword{smooth}
\references{
Silverman, B. (1986) Density Estimation, pp100-104.
}
