\name{latex.table.rq}
\alias{latex.table.rq}
\title{ Function to make a latex table from table.rq output }
\description{
	Makes a latex formatted table from the output of table.rq
	As currently configured the table has a row  for each parameter
	of the model and columns for each quantile estimated.
}
\usage{
latex.table.rq(object, transpose = FALSE, caption="caption goes here.", 
	digits=3, file=as.character(substitute(object)),\dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ table.rq \code{object}  }
  \item{transpose}{  transpose the table so quantiles are the  row dimension  }
  \item{caption}{ latex \code{caption}  }
  \item{digits}{ number of significant \code{digits} }
  \item{file}{ name of output \code{file} }
  \item{\dots}{ other optional arguments }
}
\details{
	calls latex.table() a modified version of Frank Harrell's Splus function
}
\value{
	writes a file and returns invisably.
}
\author{ Roger Koenker }


\seealso{ \code{\link{table.rq}}}

\examples{
}
\keyword{ utilities }% at least one, from doc/KEYWORDS
