\name{rq.process.object}
\alias{rq.process.object}
\title{
Linear Quantile Regression Process Object 
}
\description{
These are objects of class \code{rq.process.}
They represent the fit of a linear conditional quantile function model. 
}
\section{Generation}{
This class of objects is returned from the \code{rq}
 function 
to represent a fitted linear quantile regression model. 
}
\section{Methods}{
The \code{"rq.process"} class of objects has 
methods for the following generic 
functions: 
\code{effects}, \code{formula}
, \code{labels}
, \code{model.frame}
, \code{model.matrix}
, \code{plot}
, \code{predict}
, \code{print}
, \code{print.summary}
, \code{summary}
}
\section{Structure}{
  The following components must be included in a legitimate \code{rq.process}
  object. 

  \describe{
    \item{\code{sol}}{
      The primal solution array.  This is a (p+3)  by  J  matrix  whose  
      first  row  contains  the 'breakpoints'
      \eqn{tau_1, tau_2, \dots, tau_J},   
      of  the  quantile function, i.e. the values in [0,1] at which  the  
      solution changes,  row  two  contains  the  corresponding quantiles 
      evaluated at the mean design point, i.e. the inner product of  
      xbar  and  \eqn{b(tau_i)}, the third row contains the value of the objective
      function evaluated at the corresponding \eqn{tau_j}, and the last p rows 
      of the matrix give \eqn{b(tau_i)}.  The solution \eqn{b(tau_i)} prevails from  
      \eqn{tau_i} to \eqn{tau_i+1}.  Portnoy (1991) shows that
      \eqn{J=O_p(n \log n)}{J=O_p(n log n)}.
    }
    \item{\code{dsol}}{
      The dual solution array.  This is a
      n by J matrix  containing the  dual  solution  corresponding to sol,
      the ij-th entry is 1 if \eqn{y_i > x_i b(tau_j)}, is 0 if \eqn{y_i <  x_i
      b(tau_j)},  and is between 0 and 1 otherwise, i.e. if the
      residual is zero. See  Gutenbrunner and Jureckova(1991)
      for a detailed discussion of the statistical
      interpretation of dsol.  The use of dsol in inference is described
      in Gutenbrunner, Jureckova, Koenker, and Portnoy (1994).
    }
  }
}
\details{
  These arrays are computed by parametric linear programming methods
  using using the exterior point (simplex-type) methods of the 
  Koenker--d'Orey algorithm based on Barrodale and Roberts median
  regression algorithm.
}
\references{
  [1] Koenker, R. W. and Bassett, G. W. (1978). Regression quantiles,
  \emph{Econometrica}, \bold{46}, 33--50.

  [2] Koenker, R. W. and d'Orey (1987, 1994).
  Computing Regression Quantiles.
  \emph{Applied Statistics}, \bold{36}, 383--393, and \bold{43}, 410--414.

  [3] Gutenbrunner, C. Jureckova, J. (1991).
  Regression quantile and regression rank score process in the
  linear model and derived statistics, \emph{Annals of Statistics},
  \bold{20}, 305--330.

  [4] Gutenbrunner, C., Jureckova, J., Koenker, R. and
  Portnoy, S. (1994)  Tests of linear hypotheses based on regression  
  rank scores.   \emph{Journal of Nonparametric Statistics}, 
  (2), 307--331.

  [5]  Portnoy, S. (1991).  Asymptotic behavior of the number of regression
  quantile breakpoints, \emph{SIAM Journal of Scientific
    and  Statistical Computing}, \bold{12}, 867--883.
}
\seealso{
\code{\link{rq}}.
}
\keyword{regression}
