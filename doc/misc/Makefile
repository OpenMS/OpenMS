# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#


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
