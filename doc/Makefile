# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

#######################################################################
# HELP TARGET
help:
	@echo ""
	@echo "Main targets are:"
	@echo "  - 'make [help]'    shows this help"
	@echo "  - 'make doc'       creates the HTML documentation"
	@echo "  - 'make tutorial'  creates the tutorial as a PDF"
	@echo ""

#######################################################################
# creates the documentation
doc:
	@rm -rf doxygen/xml_output doxygen/html_output
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile

#######################################################################
# creates the internal documentation
idoc:
	@rm -rf doxygen/xml_output doxygen/html_output
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile_internal

#######################################################################
# creates the tutorial
tutorial: tutorial/Tutorial.doxygen tutorial/Doxyfile
	@rm -rf tutorial/latex_output
	@doxygen tutorial/Doxyfile
	@cp tutorial/refman_overwrite.tex tutorial/latex_output/refman.tex
	@cd tutorial/latex_output/ && make
	@cp tutorial/latex_output/refman.pdf tutorial/Tutorial.pdf
