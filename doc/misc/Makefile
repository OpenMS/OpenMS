# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# $Id: Makefile,v 1.1.1.1 2005/04/04 13:03:05 marc_sturm Exp $
# $Author: marc_sturm $


SOURCES = coding_style

PDF=$(SOURCES:=.pdf)
PS=$(SOURCES:=.ps)
DVI=$(SOURCES:=.dvi)

default: $(PDF)

all: $(PDF) $(PS) $(DVI) rights

ps: $(PS)

%.pdf: %.dvi
	dvipdf $<

%.ps: %.dvi
	dvips $< -o $*.ps

%.dvi: %.tex
	latex $<

realclean: clean
	rm -f *.ps *.pdf

clean:
	rm -f *.aux *.log *.dvi *.pdf *~
