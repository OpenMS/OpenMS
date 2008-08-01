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

include ../config_defs.mak

DIRECTORY=APPLICATIONS

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)


CPP_SOURCES=\
	TOPPBase.C\
	TOPPViewBase.C\
	INIFileEditorWindow.C

MOC_SOURCES=\
	TOPPViewBase_moc.C\
	INIFileEditorWindow_moc.C

SOURCES		= $(CPP_SOURCES) $(MOC_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)


default:	echodir $(OBJECTS)

-include .Dependencies

depend:	$(MOC_SOURCES)
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX_MAKEDEPEND) $(MAKEDEP_OPTS) $(CPP_SOURCES) $(MOC_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies

default_clean:
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES)

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/libOpenMS.objects ; done

include ../config.mak

TOPPViewBase_moc.C:  $(HEADER_DIR)/TOPPViewBase.h
	$(MOC) $(HEADER_DIR)/TOPPViewBase.h -o TOPPViewBase_moc.C
	
INIFileEditorWindow_moc.C:  $(HEADER_DIR)/INIFileEditorWindow.h
	$(MOC) $(HEADER_DIR)/INIFileEditorWindow.h -o INIFileEditorWindow_moc.C
	

# For debugging the TOPP applications, please create the file 'Makefile_debug'
# and add the line 'ADD_CXXFLAGS=-DDEBUG_TOPP' to it
-include Makefile_debug

