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

include ../config_defs.mak

DIRECTORY=APPLICATIONS

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

THISLIB=$(LIBNAME)

CPP_SOURCES=\
	TOPPBase.C

# These depend on QT:
GUI_CPP_SOURCES=\
	TOPPViewBase.C		

GUI_MOC_SOURCES=\
	TOPPViewBase_moc.C

#conditional compiles QT
ifeq ($(ENABLE_GUI),true)
	ADD_CPP_SOURCES=\
		$(GUI_CPP_SOURCES)
	ADD_MOC_SOURCES=\
		$(GUI_MOC_SOURCES)
endif

ADD_INCLUDES=$(GUI_INCLUDES)


SOURCES		= $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(ADD_MOC_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)


default:	echodir $(OBJECTS)

-include .Dependencies

depend:	$(MOC_SOURCES) $(ADD_MOC_SOURCES)
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies

default_clean:
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES) $(ADD_MOC_SOURCES)

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/$(THISLIB).objects ; done

include ../config.mak

TOPPViewBase_moc.C:  $(HEADER_DIR)/TOPPViewBase.h
	$(MOC) $(HEADER_DIR)/TOPPViewBase.h -o TOPPViewBase_moc.C

# For debugging the TOPP applications, please create the file 'Makefile_debug'
# and add the line 'ADD_CXXFLAGS=-DDEBUG_TOPP' to it
-include Makefile_debug