\name{barro}
\alias{barro}
\title{Barro Data}
\description{
Version of the Barro Growth Data used in Koenker and Machado(1999).
This is a regression data set consisting of 161 observations on determinants
of cross country GDP growth rates. There are 13 covariates with dimnames
corresponding to the original Barro and Lee source.  See
http://www.nber.org/pub/barro.lee/.  The first 71 observations are on
the period 1965-75, remainder on 1987-85.  
}
\usage{data(barro)}
\format{A data frame containing 161 observations on 14 variables:
\tabular{rl}{
        [,1] \tab "Annual Change Per Capita GDP"\cr
        [,2] \tab "Initial Per Capita GDP"\cr
        [,3] \tab "Male Secondary Education"\cr
        [,4] \tab "Female Secondary Education"\cr
        [,5] \tab "Female Higher Education"\cr
        [,6] \tab "Male Higher Education"\cr
        [,7] \tab "Life Expectancy"\cr
        [,8] \tab "Human Capital"\cr
        [,9] \tab "Education/GDP"\cr
        [,10] \tab "Investment/GDP"\cr
        [,11] \tab "Public Consumption/GDP"\cr
        [,12] \tab "Black Market Premium"\cr
        [,13] \tab "Political Instability"\cr
        [,14] \tab "Growth Rate Terms Trade"}

}
\references{
Koenker, R. and J.A.F. Machado (1999) Goodness of Fit and Related Inference Processes for Quantile Regression, JASA, 1296-1310.}
\keyword{datasets}
