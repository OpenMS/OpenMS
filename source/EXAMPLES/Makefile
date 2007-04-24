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

CPP_SOURCES=\
  Tutorial_DPosition.C\
  Tutorial_DRange.C\
	Tutorial_MSSpectrum.C\
	Tutorial_MSExperiment.C\
	Tutorial_FeatureMap.C\
	Tutorial_RangeManager.C\
	Tutorial_MetaInfo.C\
	Tutorial_FileIO.C\
	Tutorial_Param.C\
	Tutorial_PeakPickerCWT.C
	

OBJECTS=$(SOURCES:.C=.o)

all:	build_all

include ../common.mak
-include .Dependencies

EXECUTABLES=$(SOURCES:.C=)

build_all:	$(EXECUTABLES)

clean:	default_clean

dummy:

Tutorial_DPosition: Tutorial_DPosition.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_DRange: Tutorial_DRange.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_MSSpectrum: Tutorial_MSSpectrum.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_MSExperiment: Tutorial_MSExperiment.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_FeatureMap: Tutorial_FeatureMap.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_RangeManager: Tutorial_RangeManager.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_MetaInfo: Tutorial_MetaInfo.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_FileIO: Tutorial_FileIO.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_Param: Tutorial_Param.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)

Tutorial_PeakPickerCWT: Tutorial_PeakPickerCWT.o
	$(CXX) $(LDFLAGS) -o $@ $< $(GUI_LIBS)
