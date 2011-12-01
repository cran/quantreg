\name{residuals.nlrq}
\alias{residuals.nlrq}
\title{ Return residuals of an nlrq object }
\description{
Set algorithmic parameters for nlrq (nonlinear quantile regression function)
}
\usage{
\method{residuals}{nlrq}(object, type = c("response", "rho"), ...) 
}
\arguments{
  \item{object}{an `nlrq' object as returned by function `nlrq'}
  \item{type}{the type of residuals to return: "response" is the distance
  between observed and predicted values; "rho" is the weighted distance used
  to calculate the objective function in the minimisation algorithm as 
  tau * pmax(resid, 0) + (1 - tau) * pmin(resid, 0), where resid are the
  simple residuals as above (with type="response").
  }
  \item{...}{further arguments passed to or from other methods.}
}
\seealso{  \code{\link{nlrq}} }

\keyword{models}  
\keyword{regression}
\keyword{nonlinear}
