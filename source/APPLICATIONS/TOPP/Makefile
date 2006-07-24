# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#		   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

DIRECTORY=APPLICATIONS/TOPP

ifeq ($(GSL_SUPPORT),true)
GSL_EXEC=\
	MapStatistics\
	AdditiveSeries\
	MapMatcher\
#	NoiseFilter\

endif

ifeq ($(LIBSVM_SUPPORT),true)
LIBSVM_EXEC=\
	RTModel\
	RTPredict
endif

EXECUTABLES=\
	FileConverter\
	FileInfo\
	FileFilter\
	FileMerger\
	TOPPView\
	FeatureFinder\
	$(GSL_EXEC)\
	BaselineFilter\
	IDFilter\
	MascotAdapter\
	InspectAdapter\
	PeakPicker\
	UnlabeledMatcher\
	LabeledMatcher\
	Dewarper\
	$(LIBSVM_EXEC)

#	FeaturePairSplitter\
#	RawDataResampler\

#	DBExporter\
#	DBImporter\


default: all

all: $(EXECUTABLES)

# replace library-specific flags by application-
# specific compiler flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 

include ../../common.mak
-include .Dependencies

clean: default_clean
	@-$(RM) -rf .Dependencies html

ADD_CXXFLAGS= -DQT_THREAD_SUPPORT -pg -DQT_SHARED -D_REENTRANT -pipe -Wall -W -Wno-deprecated -DQT_NO_DEBUG 

OPENMS_INCLUDES=$(GUI_INCLUDES) -I./
LIBS=$(GUI_LIBS)

CPP_SOURCES =\
	TOPPBase.C\
	$(EXECUTABLES:=.C)

MapStatistics: MapStatistics.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

MapMatcher: MapMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

NoiseFilter: NoiseFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

BaselineFilter: BaselineFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

AdditiveSeries: AdditiveSeries.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

RTModel: RTModel.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

RTPredict: RTPredict.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

FileConverter: FileConverter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

FileInfo: FileInfo.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

FileFilter: FileFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

TOPPView: TOPPView.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

FeatureFinder: FeatureFinder.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

SpectrumComparator: SpectrumComparator.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

SpectrumLibBuilder: SpectrumLibBuilder.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

SpectrumLibMatcher: SpectrumLibMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

Dewarper: Dewarper.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

PeakPicker: PeakPicker.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

IDFilter: IDFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

FeaturePairSplitter: FeaturePairSplitter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

DBExporter: DBExporter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

DBImporter: DBImporter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

LabeledMatcher: LabeledMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

MascotAdapter: MascotAdapter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

InspectAdapter: InspectAdapter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

RawDataResampler: RawDataResampler.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

UnlabeledMatcher: UnlabeledMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)

FileMerger: FileMerger.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
