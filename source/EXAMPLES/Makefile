# -*- Mode: make; tab-width: 2; -*-
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
# $Id: Makefile,v 1.19 2006/06/01 16:51:16 groepl Exp $
# $Author: groepl $
# $Maintainer: Ole Schulz-Trieglaff $
# --------------------------------------------------------------------------

DIRECTORY=EXAMPLES

CPP_SOURCES= \
	FeatureFinder.C \
	DimensionDescription_example.C\
	FeatureMapping.C \
	Serialization_example_GPS.C \
	Serialization_example_KERNEL_1.C \
	Serialization_example_KERNEL_2.C

OBJECTS=$(SOURCES:.C=.o)

all:	build_all

include ../common.mak
-include .Dependencies

EXECUTABLES=$(SOURCES:.C=)

build_all:	$(EXECUTABLES)

# remove the library specific flags, includes, libs
OPENMS_INCLUDES=$(GUI_INCLUDES)
LIBS=$(GUI_LIBS)
LIB_CXXFLAGS = $(NONLIB_CXXFLAGS)

CXXFLAGS += -g

clean:	default_clean

DimensionDescription_example: DimensionDescription_example.o
	$(CXX) $(LDFLAGS) -o $@ $@.o $(LIBS)

FeatureFinder:	FeatureFinder.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

FeatureMapping: FeatureMapping.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

Serialization_example_GPS.o: Serialization_example_GPS.C

Serialization_example_GPS: Serialization_example_GPS.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

Serialization_example_KERNEL_1.o: Serialization_example_KERNEL_1.C

Serialization_example_KERNEL_1: Serialization_example_KERNEL_1.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

Serialization_example_KERNEL_2.o: Serialization_example_KERNEL_2.C

Serialization_example_KERNEL_2: Serialization_example_KERNEL_2.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

dummy:

