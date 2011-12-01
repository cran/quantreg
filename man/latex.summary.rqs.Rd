\name{latex.summary.rqs}
\alias{latex.summary.rqs}
\title{ Make a latex table from a table of rq results}
\description{
  Produces a file with latex commands for a table of rq results.
}
\usage{\method{latex}{summary.rqs}(x, transpose = FALSE, caption = "caption goes here.", 
		digits = 3, file = as.character(substitute(x)), ...)
}
\arguments{
  \item{x}{\code{x} is an object of class \code{summary.rqs}}
  \item{transpose}{if \code{TRUE} transpose table so that 
        rows are quantiles and columns are covariates. }
  \item{caption}{ caption for the table}
  \item{digits}{ decimal precision of table entries.}
  \item{file}{ name of file }
  \item{\dots}{ optional arguments for \code{latex.table}  }
}
\details{
  Calls  \code{latex.table}. 
}
\value{
  Returns invisibly after writing the file.
}
\author{ R. Koenker}
\seealso{ \code{\link{summary.rqs}}, \code{\link{latex.table}}}
\keyword{IO}
