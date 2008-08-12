# -*- mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#		   OpenMS Mass Spectrometry Framework
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

include ../../config_defs.mak

DIRECTORY=APPLICATIONS/UTILS

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

UTILS_EXECUTABLES=\
	FuzzyDiff\
	Digestor\
	IdXMLInfo\
	XMLValidator\
	FFEval

all: $(UTILS_EXECUTABLES)

# replace library-specific flags by application-specific compiler flags 
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 

include ../../common.mak
-include .Dependencies

clean: default_clean
	@rm -rf $(UTILS_EXECUTABLES)
	@-$(RM) -rf .Dependencies html

ADD_CXXFLAGS= -D_REENTRANT

CPP_SOURCES =\
	$(UTILS_EXECUTABLES:=.C)

%: %.o
	$(ECHO) $(CXX) $(LDFLAGS) -o $@ $< $(LIBS); \
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS); \
	$(ECHO) $(CP) $@ $(PREFIX)/bin/; \
	$(CP) $@ $(PREFIX)/bin/; \
