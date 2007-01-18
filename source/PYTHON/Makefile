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
# $Maintainer: $
# --------------------------------------------------------------------------
#

DIRECTORY=PYTHON

# allow processing of SIP files
.SUFFIXES: .sip

# default target: the OpenMS Python extensions in OpenMSmodule.so
all:	OpenMSmodule.so

#
include ../common.mak
-include OpenMS.sbf

SIP_SOURCES=\
	$(sources)

SIP_OBJECTS	=	\
	$(sources:.C=.o)

SIP_HEADERS = \
	$(PYOPENMS_SOURCES:.C=.h) \

SIP_HEADER_FILES =\
	$(headers)

SIP_OPTS=-c SIP -e -s ".C" -b OpenMS.sbf
SIPFLAGS += $(PYTHON_INCLUDES) $(SIP_INCLUDES)
CXXFLAGS += $(SIPFLAGS)
 
OPENMS_INCLUDES+=$(VIEW_INCLUDES)
OPENMS_LIBS=$(GUI_LIBS)
ADD_INCLUDES=$(PYTHON_INCLUDES) $(SIP_INCLUDES)
 
PYOPENMS_SOURCES = \

PYOPENMS_OBJECTS  = $(PYOPENMS_SOURCES:.C=.o)
 
CPP_SOURCES = $(SIP_SOURCES:.h=.C) $(PYOPENMS_SOURCES)
CPP_OBJECTS = $(SIP_OBJECTS) $(PYOPENMS_OBJECTS)

OpenMSmodule.so: sipOpenMScmodule.o

sip:
	@if test ! -d SIP ; then mkdir SIP 2>/dev/null ; fi
	@echo "running SIP..."
	@$(SIP) $(SIP_OPTS) OpenMS.sip
	@echo "removing SIP headers..."
	@-rm -f SIP/tmp 2>/dev/null
	@for i in SIP/*; do mv $$i SIP/tmp ; cat SIP/tmp |sed 4d > $$i ; rm -rf SIP/tmp 2>/dev/null ; done
	@echo "updating modified files..."
	@for i in SIP/*; do if cmp -s $$i `basename $$i`; then true ; else echo " - $$i"; cp $$i . ; fi ; done


OpenMSmodule.so:  $(CPP_OBJECTS)
	$(DYNAR) $(DYNAROPTS) OpenMSmodule.so $(SIP_OBJECTS) $(PYOPENMS_OBJECTS) $(SIP_LIB) -L$(OPENMS_PATH)/lib/$(BINFMT) $(OPENMS_LIBS)

install:	OpenMSmodule.so
	cp OpenMSmodule.so  $(OPENMS_PATH)/lib/$(BINFMT)
	 
CLEAN_FILES = \
	OpenMSmodule.so\
	$(CPP_OBJECTS) \
	$(SIP_SOURCES) \
	$(SIP_HEADERS) \
	$(SIP_HEADER_FILES) \
	$(SIP_OBJECTS)
 
clean:
	@echo "removing object files, libraries, and source/header files produced by SIP..."
	@-$(RM) -f $(CLEAN_FILES) 2>/dev/null
	@-$(RM) -fr SIP

depend: sip  Makefile.inc .Dependencies

-include .Dependencies
