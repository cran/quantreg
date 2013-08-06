\name{latex.table}
\alias{latex.table}
\title{ Writes a latex formatted table to a file}
\description{
	Automatically generates a latex formatted table from the matrix x
	Controls rounding, alignment, etc, etc
}
\usage{\method{latex}{table}(x, file=as.character(substitute(x)), 
	rowlabel=file, rowlabel.just="l", cgroup, n.cgroup, rgroup, n.rgroup=NULL, 
	digits, dec, rdec, cdec, append=FALSE, dcolumn=FALSE, cdot=FALSE, 
	longtable=FALSE, table.env=TRUE, lines.page=40, caption, caption.lot, 
	label=file, double.slash=FALSE,\dots)
}
\arguments{
  \item{x}{ A matrix \code{x} with dimnames}
  \item{file}{ Name of output \code{file} (.tex will be added) }
  \item{rowlabel}{ If `x' has row dimnames, rowlabel is a character 
	string containing the column heading for the row dimnames. 
	The default is the name of the argument for x. }
  \item{rowlabel.just}{ If `x' has row dimnames, specifies the 
	justification for printing them.  Possible values are `
	"l", "r", "c"'. The heading (`rowlabel') itself
	is left justified if `rowlabel.just="l"', otherwise it is centered. }
  \item{cgroup}{ a vector of character strings defining major column headings. 
	The default is to have none. }
  \item{n.cgroup}{ a vector containing the number of columns for which each 
	element in cgroup is a heading.  For example, specify `cgroup=
	c("Major 1","Major 2")', `n.cgroup=c(3,3)' if "Major 1" is to 
	span columns 1-3 and "Major 2" is to span columns 4-6.  
	`rowlabel' does not count in the column numbers.  You can omit 
	`n.cgroup' if all groups have the same number of columns.  }
  \item{rgroup}{ a vector of character strings containing headings for row 
	groups.  `n.rgroup' must be present when `rgroup' is given. The 
	first `n.rgroup[1]' rows are sectioned off and `rgroup[1]' is 
	used as a bold heading for them. The usual row dimnames (which must 
	be present if `rgroup' is) are indented. The next `n.rgroup[2]' 
	rows are treated likewise, etc.  }
  \item{n.rgroup}{integer vector giving the number of rows in each grouping. 
	If `rgroup' is not specified, `n.rgroup' is just used to divide off 
	blocks of rows by horizontal lines. If `rgroup' is given but 
	`n.rgroup' is omitted, `n.rgroup' will default so that each row 
	group contains the same number of rows.  }
  \item{digits}{ causes all values in the table to be formatted to `digits' 
	significant digits.  `dec' is usually preferred.  }
  \item{dec}{ If `dec' is a scalar, all elements of the matrix will be 
	rounded to `dec' decimal places to the right of the decimal. 
	`dec' can also be a matrix whose elements correspond to `x', for 
	customized rounding of each element.}
  \item{rdec}{ a vector specifying the number of decimal places to the right 
	for each row (`cdec' is more commonly used than `rdec') }
  \item{cdec}{ a vector specifying the number of decimal places for each 
	column }
  \item{append}{  defaults to `F'. Set to `T' to append output to an 
	existing file.}
  \item{dcolumn}{ Set to `T' to use David Carlisles `dcolumn' style for 
	decimal alignment.
	Default is `F', which aligns columns of numbers by changing leading
	blanks to "~", the LaTeX space-holder. You will probably want to
	use `dcolumn' if you use `rdec', as a column may then contain varying
	number of places to the right of the decimal. `dcolumn' can line up
	all such numbers on the decimal point, with integer values right-
	justified at the decimal point location of numbers that actually
	contain decimal places.  }
  \item{cdot}{ Set to `T' to use centered dots rather than ordinary periods 
	in numbers.}
  \item{longtable}{ Set to `T' to use David Carlisles LaTeX `longtable' style, 
	allowing long tables to be split over multiple pages with headers 
	repeated on each page.}
  \item{table.env}{Set `table.env=FALSE' to suppress enclosing the table
	in a LaTeX `table' environment.  `table.env' only applies when 
	`longtable=FALSE'.  You may not specify a `caption' if `table.env=FALSE'.  }
  \item{lines.page}{ Applies if `longtable=TRUE'. No more than `lines.page' 
	lines in the body of a table will be placed on a single page. 
	Page breaks will only occur at `rgroup' boundaries. }
  \item{caption}{ a text string to use as a caption to print at the top of the 
	first page of the table. Default is no caption.  }
  \item{caption.lot}{ a text string representing a short caption to be used 
	in the "List of Tables".  By default, LaTeX will use `caption'.  }
  \item{label}{ a text string representing a symbolic label for the table 
	for referencing with the LaTex `\\ref\{label\}' command. The default 
	is `file'.  `label' is only used if `caption' is given.  }
  \item{double.slash}{set to `T' to output `\\' as `\\\\' in LaTeX commands. 
	Useful when you are reading the output file back into an S vector 
	for later output.  }
  \item{\dots}{other optional arguments}
}
\value{ returns invisibly }
\references{ Minor modification of Frank Harrell's Splus code }
\author{ Roger Koenker }
\keyword{ utilities }
