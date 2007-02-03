# -*- Mode: make; tab-width: 2; -*-
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
# $Maintainer: Ole Schulz-Trieglaff $
# --------------------------------------------------------------------------

DIRECTORY=EXAMPLES

CPP_SOURCES= \
	Serialization_example_GPS.C \
	Serialization_example_KERNEL_1.C \
	FeatureMapping.C \

### Apparently this one doesn't compile with g++ version 4.1.0.
### It used to compile with an earlier version...
### TODO fix that (Clemens).
#	Serialization_example_KERNEL_2.C \


OBJECTS=$(SOURCES:.C=.o)

all:	build_all

include ../common.mak
-include .Dependencies

EXECUTABLES=$(SOURCES:.C=)

build_all:	$(EXECUTABLES)

OPENMS_INCLUDES =  $(GUI_INCLUDES) 

LIBS = $(GUI_LIBS)

LIB_CXXFLAGS = $(NONLIB_CXXFLAGS)

CXXFLAGS += -g

clean:	default_clean

FeatureFinder:	FeatureFinder.o 
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

FeatureMapping: FeatureMapping.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)


# I admit that this is a weird hack, but we need to prefer "boost-1_33_1/boost" over "boost".
# The .C.o rule is otherwise the same as in config.mak, which turns out not to be affected by configure (momentarily).
# (Maybe there is a problem with the contrib package?)
Serialization_example_GPS.o: Serialization_example_GPS.C
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(LIB_CXXFLAGS) $(OPENMS_INCLUDES:include=include/boost-1_33_1) $(OPENMS_INCLUDES) $(ADD_INCLUDES) $(GUI_INCLUDES) -c $*.C -o $@

Serialization_example_GPS: Serialization_example_GPS.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)


# I admit that this is a weird hack, but we need to prefer "boost-1_33_1/boost" over "boost".
# The .C.o rule is otherwise the same as in config.mak, which turns out not to be affected by configure (momentarily).
# (Maybe there is a problem with the contrib package?)
Serialization_example_KERNEL_1.o: Serialization_example_KERNEL_1.C
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(LIB_CXXFLAGS) $(OPENMS_INCLUDES:include=include/boost-1_33_1) $(OPENMS_INCLUDES) $(ADD_INCLUDES) $(GUI_INCLUDES) -c $*.C -o $@

Serialization_example_KERNEL_1: Serialization_example_KERNEL_1.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)


# I admit that this is a weird hack, but we need to prefer "boost-1_33_1/boost" over "boost".
# The .C.o rule is otherwise the same as in config.mak, which turns out not to be affected by configure (momentarily).
# (Maybe there is a problem with the contrib package?)
Serialization_example_KERNEL_2.o: Serialization_example_KERNEL_2.C
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(LIB_CXXFLAGS) $(OPENMS_INCLUDES:include=include/boost-1_33_1) $(OPENMS_INCLUDES) $(ADD_INCLUDES) $(GUI_INCLUDES) -c $*.C -o $@

Serialization_example_KERNEL_2: Serialization_example_KERNEL_2.o
	$(CXX)	$(LDFLAGS) -o $@ $@.o $(LIBS)

dummy:


bla:
	echo $(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(LIB_CXXFLAGS) $(OPENMS_INCLUDES) $(ADD_INCLUDES) $(GUI_INCLUDES)
