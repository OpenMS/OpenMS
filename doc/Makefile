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

-include ../source/config.mak

#######################################################################
# HELP TARGET
help:
	@echo ""
	@echo "Main targets are:"
	@echo "  - 'make [help]'           shows this help"
	@echo "  - 'make doc'              creates the HTML documentation (including HTML tutorials)"
	@echo "  - 'make OpenMS_tutorial'  creates the OpenMS tutorial as a PDF"
	@echo "  - 'make TOPP_tutorial'    creates the TOPP tutorial as a PDF"
	@echo "  - 'make clean'            removes doxygen docu and tutorial"
	@echo ""

#######################################################################
# create parameter documentation of classes derived from DefaultParamHandler 
paramdoc: doxygen/parameters/DefaultParamHandlerDocumenter.C ../source/config.mak dummy
	@echo "";
	@echo "Note: In order to build the parameter documentation for most of the algorithms, a functioning OpenMS installation and a running X-server is required!";
	@echo "";
	@cd doxygen/parameters/ && \
	 $(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(LIB_CXXFLAGS) $(OPENMS_INCLUDES) $(ADD_INCLUDES) -c DefaultParamHandlerDocumenter.C -o DefaultParamHandlerDocumenter.o && \
	 $(CXX) $(LDFLAGS) -o DefaultParamHandlerDocumenter DefaultParamHandlerDocumenter.o $(GUI_LIBS) && \
	 ./DefaultParamHandlerDocumenter;
	@touch doxygen/parameters/DefaultParameters.doxygen

#######################################################################
# creates the documentation
doc:
	@make -i paramdoc
	@rm -rf doxygen/xml_output doxygen/html_output
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile

#######################################################################
# creates the internal documentation
idoc:
	@make -i paramdoc
	@rm -rf doxygen/xml_output doxygen/html_output
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile_internal

#######################################################################
# creates the OpenMS tutorial
OpenMS_tutorial: OpenMS_tutorial/OpenMS_Tutorial.doxygen OpenMS_tutorial/Doxyfile dummy
	@rm -rf OpenMS_tutorial/latex_output
	@doxygen OpenMS_tutorial/Doxyfile
	@cp OpenMS_tutorial/refman_overwrite.tex OpenMS_tutorial/latex_output/refman.tex
	@cd OpenMS_tutorial/latex_output/ && make
	@cp OpenMS_tutorial/latex_output/refman.pdf OpenMS_tutorial.pdf


#######################################################################
# creates the TOPP tutorial
TOPP_tutorial: TOPP_tutorial/TOPP_Tutorial.doxygen TOPP_tutorial/Doxyfile dummy
	@rm -rf TOPP_tutorial/latex_output
	@doxygen TOPP_tutorial/Doxyfile
	@cp TOPP_tutorial/refman_overwrite.tex TOPP_tutorial/latex_output/refman.tex
	@cd TOPP_tutorial/latex_output/ && make
	@cp TOPP_tutorial/latex_output/refman.pdf TOPP_tutorial.pdf

#######################################################################
# cleaning up
clean:
	@rm -rf doxygen/html_output/ doxygen/xml_output/ doxygen/doxygen-error.log
	@rm -rf OpenMS_tutorial/latex_output/ OpenMS_tutorial.pdf OpenMS_tutorial/doxygen-error.log
	@rm -rf TOPP_tutorial/latex_output/ TOPP_tutorial.pdf TOPP_tutorial/doxygen-error.log
	@rm -rf doxygen/parameters/DefaultParamHandlerDocumenter doxygen/parameters/DefaultParamHandlerDocumenter.o

dummy:
