# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#		   OpenMS Mass Spectrometry Framework
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

include ../../config_defs.mak

DIRECTORY=APPLICATIONS/TOPP

ifeq ($(CGAL_SUPPORT),true)
CGAL_EXEC=\
	UnlabeledMatcher\
	MapAlignment
endif

ifeq ($(LIBSVM_SUPPORT),true)
LIBSVM_EXEC=\
	RTModel\
	RTPredict
endif

EXECUTABLES=\
	$(CGAL_EXEC)\
	$(LIBSVM_EXEC)\
	FileConverter\
	FileInfo\
	FileFilter\
	FileMerger\
	TOPPView\
	BaselineFilter\
	IDFilter\
	MascotAdapter\
	InspectAdapter\
	SequestAdapter\
	PeakPicker\
	LabeledMatcher\
	MapDewarper\
	FeaturePairSplitter\
	DTAExtractor\
	DBExporter\
	DBImporter\
	AnalysisXMLMerger\
	MapStatistics\
	AdditiveSeries\
	MapMatcher\
	NoiseFilter\
	SpectraFilter\
	Resampler\
	FeatureFinder\
	PILISModel\
	PILISIdentification\
	OMSSAAdapter


all: $(EXECUTABLES)

# replace library-specific flags by application-
# specific compiler flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 

include ../../common.mak
-include .Dependencies

clean: default_clean
	@-$(RM) -rf .Dependencies html 
	@cd $(OPENMS_PATH)/bin/ 2>/dev/null | cat >/dev/null && rm -rf $(EXECUTABLES)

ADD_CXXFLAGS= -D_REENTRANT

OPENMS_INCLUDES=$(GUI_INCLUDES)
LIBS=$(GUI_LIBS)

CPP_SOURCES =\
	$(EXECUTABLES:=.C)

test:
	cd ${OPENMS_PATH}/source/TEST/TOPP/ && make

MapStatistics: MapStatistics.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/             

MapMatcher: MapMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

NoiseFilter: NoiseFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

BaselineFilter: BaselineFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

AdditiveSeries: AdditiveSeries.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

RTModel: RTModel.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

RTPredict: RTPredict.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

FileConverter: FileConverter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

FileInfo: FileInfo.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

FileFilter: FileFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

TOPPView: TOPPView.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

FeatureFinder: FeatureFinder.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

SpectrumComparator: SpectrumComparator.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

SpectrumLibBuilder: SpectrumLibBuilder.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

SpectrumLibMatcher: SpectrumLibMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

MapDewarper: MapDewarper.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

PeakPicker: PeakPicker.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

IDFilter: IDFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

FeaturePairSplitter: FeaturePairSplitter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

DBExporter: DBExporter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

DBImporter: DBImporter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

LabeledMatcher: LabeledMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

MascotAdapter: MascotAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

InspectAdapter: InspectAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

SequestAdapter: SequestAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/

Resampler: Resampler.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

UnlabeledMatcher: UnlabeledMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

FileMerger: FileMerger.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

SpectraFilter: SpectraFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

DTAExtractor: DTAExtractor.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/

MapAlignment: MapAlignment.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

AnalysisXMLMerger: AnalysisXMLMerger.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/   

PILISIdentification: PILISIdentification.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/

PILISModel: PILISModel.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/

OMSSAAdapter: OMSSAAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(OPENMS_PATH)/bin/

# For debugging the TOPP applications, please create the file 'Makefile_debug'
# in source/APPLICATIONS/ and add the line 'ADD_CXXFLAGS=-DDEBUG_TOPP' to it
-include ../Makefile_debug
