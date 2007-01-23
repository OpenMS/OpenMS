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
# $Maintainer:  Long John Silver$
# --------------------------------------------------------------------------

include ../../config_defs.mak

DIRECTORY=FORMAT/DB

THISLIB=$(LIBNAME)

CPP_SOURCES=

ifeq ($(ENABLE_DB),true)
	DB_CPP_SOURCES=\
		DBConnection.C\
		DBAdapter.C
endif

ADD_INCLUDES=$(GUI_INCLUDES)

SOURCES		= $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(DB_CPP_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

default:	echodir $(OBJECTS)

-include .Dependencies

depend:
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(DB_CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies

default_clean:
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS)

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/$(THISLIB).objects ; done

include ../../config.mak