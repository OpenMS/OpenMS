# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2004 -- Oliver Kohlbacher, Knut Reinert
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
# $Id: Makefile,v 1.3 2006/03/08 17:34:21 groepl Exp $
# $Author: groepl $
# $Maintainer: Clemens Groepl $
# --------------------------------------------------------------------------
#

EXECUTABLES=\
	SimpleFeatureMatcher	
# FeatureMatcher

DIRECTORY=APPLICATIONS/FEATUREMATCHER

default: all

all: $(EXECUTABLES)

# replace library-specific flags by application-
# specific compiler flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 

include ../../common.mak
-include .Dependencies

clean: default_clean
	@-$(RM) -rf .Dependencies


ADD_CXXFLAGS= -DQT_THREAD_SUPPORT -DQT_SHARED -D_REENTRANT -pipe -Wall -W -O2 -Wno-deprecated -DQT_NO_DEBUG

# suppress annoying warnings from Qt XML header files
CPP_MODE_FLAGS+= -Wno-non-virtual-dtor

OPENMS_INCLUDES=$(GUI_INCLUDES)
LIBS=$(GUI_LIBS)

CPP_SOURCES= \
	FeatureMatcher.C \
	SimpleFeatureMatcher.C 

FeatureMatcher:	FeatureMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $@.o $(LIBS)

SimpleFeatureMatcher: SimpleFeatureMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $@.o $(LIBS)

# debugging
usual:
	@echo +++++ \# make: target \"usual\" BEGIN ++++++++++++++
	rm -f TOPP.log bla.out
	SimpleFeatureMatcher -ini ../TOPP/TOPP.ini -d 1 -in1 features_04111706.xml -in2 features_04111707.xml -out bla.out
	@echo --------------- contents of TOPP.log: ---------------
	cat TOPP.log
	@echo +++++ \# make: target \"usual\" END ++++++++++++++++
