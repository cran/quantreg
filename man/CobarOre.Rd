\name{CobarOre}
\alias{CobarOre}
\docType{data}
\title{ Cobar Ore data }
\description{
Cobar Ore data from Green and Silverman (1994). 
The data consists of measurements on the "true width"
of an ore-bearing rock layer from a mine in Cobar, Australia.
}
\usage{data(CobarOre)}
\format{
  A data frame with 38 observations on the following 3 variables.
  \describe{
    \item{x}{x-coordinate of location of mine site}
    \item{y}{y-coordinate of location of mine site}
    \item{z}{ore thickness}
  }
}
\source{
Green, P.J. and B.W. Silverman (1994) Nonparametric Regression Generalized Linear Models:
A roughness penalty approach, Chapman Hall.
}
\examples{
data(CobarOre)
plot(CobarOre)
}
\keyword{datasets}
