# -*- mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
	@echo "  - 'make doc'              creates the HTML documentation containing installation"
	@echo "                            instructions, API documentation and tutorials."
	@echo "  - 'make OpenMS_tutorial'  creates the OpenMS tutorial in PDF format"
	@echo "  - 'make TOPP_tutorial'    creates the TOPP tutorial in PDF format"
	@echo "  - 'make clean'            removes doxygen docu and tutorial"
	@echo ""

#######################################################################
# create parameter documentation of classes derived from DefaultParamHandler 
paramdoc: dummy
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "Creating the algorithm parameter documentation";
	@echo "";
	@echo "Note: A functioning OpenMS installation and a running X-server is required for this step!";
	@echo "      If this step fails the rest of the documentation is created nevertheless.";
	@echo "";
	@-make paramdoc_try 

# the actual paramdoc is created in this target to make sure the above text is printed without the config.mak file
paramdoc_try: doxygen/parameters/DefaultParamHandlerDocumenter.C ../source/config.mak dummy
	@cd doxygen/parameters/ && \
	 $(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(LIB_CXXFLAGS) $(OPENMS_INCLUDES) $(ADD_INCLUDES) -c DefaultParamHandlerDocumenter.C -o DefaultParamHandlerDocumenter.o && \
	 $(CXX) $(LDFLAGS) -o DefaultParamHandlerDocumenter DefaultParamHandlerDocumenter.o $(LIBS) && \
	 ./DefaultParamHandlerDocumenter;
	@echo "[ DefaultParamHandlerDocumenter: ok ]"
	@touch doxygen/parameters/DefaultParameters.doxygen

#######################################################################
# creates the documentation
doc:
	@make paramdoc
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "Creating the documentation";
	@echo "";
	@rm -rf doxygen/xml_output html
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "The documentation has been successfully created.";
	@echo "You can now open 'index.html' in a web browser.";
	@echo "";

#######################################################################
# creates the internal documentation
idoc:
	@make paramdoc
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "Creating the internal documentation";
	@echo "";
	@rm -rf doxygen/xml_output html
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile_internal
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "The documentation has been successfully created.";
	@echo "You can now open 'index.html' in a web browser.";
	@echo "";

#######################################################################
# Creates the full-blown, graphical, super-duper, internal documentation.
# This requires the "dot" program, see http://www.graphviz.org/
# Also required is some patience when creating the doc.
# Therefore the output goes to a different directory.
dotdoc:
	@make paramdoc
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "Creating the DOT documentation";
	@echo "";
	@rm -rf doxygen/xml_output html-dot
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile_dot
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "The documentation has been successfully created.";
	@echo "You can now open 'index.html' in a web browser.";
	@echo "";


#######################################################################
# Creates the internal documentation without class documentation.
noclassdoc:
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "Creating the internal documentation with out class documentation";
	@echo "";
	@rm -rf doxygen/xml_output html
	@test -f doxygen/Doxyfile || cp ../source/config/Doxyfile.in doxygen/Doxyfile
	@doxygen doxygen/Doxyfile_noclass
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "The documentation has been successfully created.";
	@echo "You can now open 'index.html' in a web browser.";
	@echo "";

#######################################################################
# creates the OpenMS tutorial
OpenMS_tutorial: OpenMS_tutorial/OpenMS_Tutorial.doxygen OpenMS_tutorial/Doxyfile dummy
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "Creating the OpenMS tutorial";
	@echo "";
	@rm -rf OpenMS_tutorial/latex_output
	@doxygen OpenMS_tutorial/Doxyfile
	@cp OpenMS_tutorial/refman_overwrite.tex OpenMS_tutorial/latex_output/refman.tex
	@cd OpenMS_tutorial/latex_output/ && make
	@cp OpenMS_tutorial/latex_output/refman.pdf OpenMS_tutorial.pdf
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "The OpenMS tutorial in PDF format has been successfully created:";
	@echo "OpenMS_tutorial.pdf";
	@echo "";


#######################################################################
# creates the TOPP tutorial
TOPP_tutorial: TOPP_tutorial/TOPP_Tutorial.doxygen TOPP_tutorial/TOPPView_Tutorial.doxygen TOPP_tutorial/Doxyfile dummy
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "Creating the TOPP tutorial (TOPP_tutorial.pdf)";
	@echo "";
	@rm -rf TOPP_tutorial/latex_output
	@doxygen TOPP_tutorial/Doxyfile
	@cp TOPP_tutorial/refman_overwrite.tex TOPP_tutorial/latex_output/refman.tex
	@cd TOPP_tutorial/latex_output/ && make
	@cp TOPP_tutorial/latex_output/refman.pdf TOPP_tutorial.pdf
	@echo "";
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	@echo "The TOPP tutorial in PDF format has been successfully created:";
	@echo "TOPP_tutorial.pdf";
	@echo "";

#######################################################################
# cleaning up
clean:
	@rm -rf html html-dot doxygen/xml_output/ doxygen/doxygen-error.log
	@rm -rf OpenMS_tutorial/latex_output/ OpenMS_tutorial.pdf OpenMS_tutorial/doxygen-error.log
	@rm -rf TOPP_tutorial/latex_output/ TOPP_tutorial.pdf TOPP_tutorial/doxygen-error.log
	@rm -rf doxygen/parameters/DefaultParamHandlerDocumenter doxygen/parameters/DefaultParamHandlerDocumenter.o

dummy:
