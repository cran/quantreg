# Makefile for LaTeX files

LATEX	= pdflatex
BIBTEX	= bibtex

SRC	= rq.tex 

all 	: rq.pdf 

%.pdf:	%.tex
	$(LATEX) $*
	$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*
	rm -f Rplots* *.ps rq-* $*.dvi $*.log $*.bbl *.tex $*.aux $*.blg $*.log
