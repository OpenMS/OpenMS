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

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

ifeq ($(LIBSVM_SUPPORT),true)
LIBSVM_EXEC=\
	RTModel\
	RTPredict
endif

TOPP_EXECUTABLES=\
	$(LIBSVM_EXEC)\
	MapAlignment\
	ConsensusID\
	FileConverter\
	FileInfo\
	FileFilter\
	FileMerger\
	BaselineFilter\
	MascotAdapter\
	InspectAdapter\
	SequestAdapter\
	PepNovoAdapter\
	PeakPicker\
	LabeledMatcher\
	MapDewarper\
	FeaturePairSplitter\
	DTAExtractor\
	DBExporter\
	DBImporter\
	IDFilter\
	IDMerger\
	AdditiveSeries\
	MapMatcher\
	NoiseFilter\
	SpectraFilter\
	Resampler\
	FeatureFinder\
	PILISModel\
	PILISIdentification\
	OMSSAAdapter\
	UnlabeledMatcher\
	MapNormalizer\
	Decharger\
	InternalCalibration\
	ExternalCalibration\
	XTandemAdapter

OTHER_EXECUTABLES=\
	INIFileEditor\
	TOPPView

all: $(TOPP_EXECUTABLES) $(OTHER_EXECUTABLES)

write_ini:
	@mkdir -p ini
	@cd ini
	@for i in $(TOPP_EXECUTABLES); do echo $$i; `cd ini && $$i -write_ini $$i.ini`; done;

# replace library-specific flags by application-
# specific compiler flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 

include ../../common.mak
-include .Dependencies

clean: default_clean
	@rm -rf $(TOPP_EXECUTABLES) $(OTHER_EXECUTABLES)
	@-$(RM) -rf .Dependencies html 
	@cd $(PREFIX)/bin/ 2>/dev/null | cat >/dev/null && rm -rf $(TOPP_EXECUTABLES) $(OTHER_EXECUTABLES)

ADD_CXXFLAGS= -D_REENTRANT

LIBS=$(GUI_LIBS)

CPP_SOURCES =\
	$(TOPP_EXECUTABLES:=.C)\
	$(OTHER_EXECUTABLES:=.C)

test:
	cd ${OPENMS_PATH}/source/TEST/TOPP/ && make

ConsensusID: ConsensusID.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/ 

MapMatcher: MapMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

NoiseFilter: NoiseFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

BaselineFilter: BaselineFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

AdditiveSeries: AdditiveSeries.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

RTModel: RTModel.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

RTPredict: RTPredict.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

FileConverter: FileConverter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

FileInfo: FileInfo.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

FileFilter: FileFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

TOPPView: TOPPView.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

IdXMLConverter: IdXMLConverter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

FeatureFinder: FeatureFinder.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

SpectrumComparator: SpectrumComparator.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

SpectrumLibBuilder: SpectrumLibBuilder.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

SpectrumLibMatcher: SpectrumLibMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

MapDewarper: MapDewarper.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

PeakPicker: PeakPicker.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

IDFilter: IDFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

FeaturePairSplitter: FeaturePairSplitter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

DBExporter: DBExporter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

DBImporter: DBImporter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

LabeledMatcher: LabeledMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

MascotAdapter: MascotAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

InspectAdapter: InspectAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

SequestAdapter: SequestAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

PepNovoAdapter: PepNovoAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

Resampler: Resampler.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

UnlabeledMatcher: UnlabeledMatcher.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

FileMerger: FileMerger.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

SpectraFilter: SpectraFilter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

DTAExtractor: DTAExtractor.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

MapAlignment: MapAlignment.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

IDMerger: IDMerger.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/   

PILISIdentification: PILISIdentification.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

PILISModel: PILISModel.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

OMSSAAdapter: OMSSAAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

XTandemAdapter: XTandemAdapter.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

Decharger: Decharger.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

MapNormalizer: MapNormalizer.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

INIFileEditor: INIFileEditor.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

InternalCalibration: InternalCalibration.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/

ExternalCalibration: ExternalCalibration.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)
	$(CP) $@ $(PREFIX)/bin/


# For debugging the TOPP applications, please create the file 'Makefile_debug'
# in source/APPLICATIONS/ and add the line 'ADD_CXXFLAGS=-DDEBUG_TOPP' to it
-include ../Makefile_debug

