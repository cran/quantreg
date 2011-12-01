\name{rq.object}
\alias{rq.object}
\alias{formula.rq}
\alias{logLik.rq}
\alias{logLik.rqs}
\alias{AIC.rq}
\alias{AIC.rqs}
\alias{extractAIC.rq}
\title{
Linear Quantile Regression Object 
}
\description{
  These are objects of class \code{"rq"}.
  They represent the fit of a linear conditional quantile function model. 
}
\section{Generation}{
  This class of objects is returned from the \code{rq} function 
  to represent a fitted linear quantile regression model. 
}
\section{Methods}{
  The \code{"rq"} class of objects has methods for the following generic 
functions: 
\code{coef}, \code{effects}
, \code{formula}
, \code{labels}
, \code{model.frame}
, \code{model.matrix}
, \code{plot}
, \code{logLik}
, \code{AIC}
, \code{extractAIC}
, \code{predict}
, \code{print}
, \code{print.summary}
, \code{residuals}
, \code{summary}
}
\section{Structure}{
  The following components must be included in a legitimate \code{rq} object. 
  \describe{
    \item{\code{coefficients}}{
      the coefficients of the quantile regression fit. 
      The names of the coefficients are the names of the 
      single-degree-of-freedom effects (the columns of the 
      model matrix). 
      If the model was fitted by method \code{"br"} with \code{ci=TRUE}, then
      the coefficient component consists of a matrix whose
      first column consists of the vector of estimated coefficients
      and the second and third columns are the lower and upper
      limits of a confidence interval for the respective coefficients.
    }
    \item{\code{residuals}}{
      the residuals from the fit. 
    }
    \item{\code{dual}}{
      the vector dual variables from the fit. 
    }
    \item{\code{rho}}{
      The value(s) of objective function at the solution.
    }
    \item{\code{contrasts}}{
      a list containing sufficient information to construct the contrasts 
      used to fit any factors occurring in the model. 
      The list contains entries that are either matrices or character vectors. 
      When a factor is coded by contrasts, the corresponding contrast matrix 
      is stored in this list. 
      Factors that appear only as dummy variables and variables in the model 
      that are matrices correspond to character vectors in the list. 
      The character vector has the level names for a factor or the column 
      labels for a matrix. 
    }
    \item{\code{model}}{
      optionally the model frame, if \code{model=TRUE}. 
    }
    \item{\code{x}}{
      optionally the model matrix, if \code{x=TRUE}. 
    }
    \item{\code{y}}{
      optionally the response, if \code{y=TRUE}. 
    }
  }
}
\details{
  The coefficients, residuals, and effects may be extracted 
  by the generic functions of the same name, rather than 
  by the \code{$} operator.   For pure \code{rq} objects this is less critical 
  than for some of the inheritor classes.  In particular, for fitted rq objects
  using "lasso" and "scad" penalties, \code{logLik} and \code{AIC} functions
  compute degrees of freedom of the fitted model as the number of estimated
  parameters whose absolute value exceeds a threshold \code{edfThresh}.  By
  default this threshold is 0.0001, but this can be passed via the \code{AIC}
  function if this value is deemed unsatisfactory.  The function \code{AIC}
  is a generic function in R, with parameter \code{k} that controls the form
  of the penalty:  the default value of \code{k} is 2 which yields the classical
  Akaike form of the penalty, while \code{k <= 0} yields the Schwarz (BIC)
  form of the penalty.
  Note that the extractor function \code{coef} returns a vector with missing values 
  omitted.  
}
\seealso{
  \code{\link{rq}},  \code{\link{coefficients}}.   
}
\keyword{regression}
