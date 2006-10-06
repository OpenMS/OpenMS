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
	NoiseFilter\
	SpectraFilter
endif

ifeq ($(LIBSVM_SUPPORT),true)
LIBSVM_EXEC=\
	RTModel\
	RTPredict
endif

EXECUTABLES=\
	$(GSL_EXEC)\
	$(LIBSVM_EXEC)\
	FileConverter\
	FileInfo\
	FileFilter\
	FileMerger\
	TOPPView\
	FeatureFinder\
	BaselineFilter\
	IDFilter\
	MascotAdapter\
	InspectAdapter\
	PeakPicker\
	UnlabeledMatcher\
	LabeledMatcher\
	Dewarper\
	FeaturePairSplitter\
	DTAExtractor\
	DBExporter\
	DBImporter

#	RawDataResampler\

all: $(EXECUTABLES)

# replace library-specific flags by application-
# specific compiler flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 

include ../../common.mak
-include .Dependencies

clean: default_clean
	@-$(RM) -rf .Dependencies html 
	@cd $(OPENMS_PATH)/bin/ 2>/dev/null | cat >/dev/null && rm -rf $(EXECUTABLES)

ADD_CXXFLAGS= -DQT_THREAD_SUPPORT -pg -DQT_SHARED -D_REENTRANT -pipe -Wall -W -Wno-deprecated -DQT_NO_DEBUG 

OPENMS_INCLUDES=$(GUI_INCLUDES) -I./
LIBS=$(GUI_LIBS)

CPP_SOURCES =\
	TOPPBase.C\
	$(EXECUTABLES:=.C)

MapStatistics: MapStatistics.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/             

MapMatcher: MapMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

NoiseFilter: NoiseFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

BaselineFilter: BaselineFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

AdditiveSeries: AdditiveSeries.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

RTModel: RTModel.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

RTPredict: RTPredict.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

FileConverter: FileConverter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

FileInfo: FileInfo.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

FileFilter: FileFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

TOPPView: TOPPView.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

FeatureFinder: FeatureFinder.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

SpectrumComparator: SpectrumComparator.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

SpectrumLibBuilder: SpectrumLibBuilder.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

SpectrumLibMatcher: SpectrumLibMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

Dewarper: Dewarper.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

PeakPicker: PeakPicker.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

IDFilter: IDFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

FeaturePairSplitter: FeaturePairSplitter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

DBExporter: DBExporter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

DBImporter: DBImporter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

LabeledMatcher: LabeledMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

MascotAdapter: MascotAdapter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

InspectAdapter: InspectAdapter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

RawDataResampler: RawDataResampler.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

UnlabeledMatcher: UnlabeledMatcher.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

FileMerger: FileMerger.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

SpectraFilter: SpectraFilter.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   

DTAExtractor: DTAExtractor.o TOPPBase.o
	$(CXX) $(LDFLAGS) -o $@ $< TOPPBase.o $(LIBS)
	cp $@ $(OPENMS_PATH)/bin/   
