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

include ../../config_defs.mak

DIRECTORY=APPLICATIONS/TOPP

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

TOPP_EXECUTABLES=\
	RTModel\
	RTPredict\
	PTModel\
	PTPredict\
	MapAligner\
	FeatureLinker\
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
	DTAExtractor\
	DBExporter\
	DBImporter\
	IDFilter\
	IDMerger\
	AdditiveSeries\
	NoiseFilter\
	SpectraFilter\
	Resampler\
	FeatureFinder\
	OMSSAAdapter\
	MapNormalizer\
	Decharger\
	InternalCalibration\
	TOFCalibration\
	TextExporter\
	FeatureFinderMRM\
	XTandemAdapter\
	FalseDiscoveryRate\
	IDDecoyProbability\
	MascotAdapterOnline\
	SILACAnalyzer\
	ITRAQAnalyzer\
	IDMapper\
	IDRTCalibration

# Those executables not derived from TOPPBase (they are left out when writing the ini files)
OTHER_EXECUTABLES=\
	INIFileEditor\
	TOPPView

all: $(TOPP_EXECUTABLES) $(OTHER_EXECUTABLES)

write_ini:
	@mkdir -p ini
	@cd ini
	@for i in $(TOPP_EXECUTABLES); do echo $$i; `cd ini && $$i -write_ini $$i.ini`; done;


write_wsdl:
	@mkdir -p wsdl
	@cd wsdl
	@for i in $(TOPP_EXECUTABLES); do echo $$i; `cd wsdl && $$i -write_wsdl $$i.xml`; done;

# replace library-specific flags by application-specific compiler flags 
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 

include ../../common.mak
-include .Dependencies

clean: default_clean
	@rm -rf $(TOPP_EXECUTABLES) $(OTHER_EXECUTABLES)
	@-$(RM) -rf .Dependencies html 
	@cd $(PREFIX)/bin/ 2>/dev/null | cat >/dev/null && rm -rf $(TOPP_EXECUTABLES) $(OTHER_EXECUTABLES)

ADD_CXXFLAGS= -D_REENTRANT

CPP_SOURCES =\
	$(TOPP_EXECUTABLES:=.C)\
	$(OTHER_EXECUTABLES:=.C)

test:
	cd ${OPENMS_PATH}/source/TEST/TOPP/ && make

splash.C:
	convert TOPPView_Splashscreen.png xpm:splash.C && mv splash.C tmp.xpm && cat tmp.xpm | sed 's/static char/const char/' > splash.C && rm tmp.xpm

### for all the stuff below, Chris is responsible. If you want to have a look at the old Makefile, see SVN-Revision 2501

### this target is only called for Windows: create an object file, from an icon resource file
%_icon.o:
	windres $@.rc $@

### this makes all the TOPP tools, TOPPView and InifileEditor and optionally includes the icon object on WindowsOS'

# special target for GUI apps, as we want to hide the additional & ugly cmd.exe screen on Windows

TOPPView: TOPPView.o splash.o
	@if [ -f $@_icon.ico -a "$(OS)" = "MINGW32" ]; then \
		$(ECHO) "Windows & Icon detected ... embedding icon into binary."; \
		$(MAKE) $@_icon.o; \
		$(ECHO) $(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $^ $@_icon.o $(LIBS); \
		$(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $^ $@_icon.o $(LIBS); \
		$(ECHO) $(CP) $@ $(PREFIX)/bin/; \
		$(CP) $@ $(PREFIX)/bin/; \
	else \
		$(ECHO) $(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $^ $(LIBS); \
		$(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $^ $(LIBS); \
		$(ECHO) $(CP) $@ $(PREFIX)/bin/; \
		$(CP) $@ $(PREFIX)/bin/; \
	fi;

INIFileEditor: INIFileEditor.o
	@if [ -f $@_icon.ico -a "$(OS)" = "MINGW32" ]; then \
		$(ECHO) "Windows & Icon detected ... embedding icon into binary."; \
		$(MAKE) $@_icon.o; \
		$(ECHO) $(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $@.o $@_icon.o $(LIBS); \
		$(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $@.o $@_icon.o $(LIBS); \
		$(ECHO) $(CP) $@ $(PREFIX)/bin/; \
		$(CP) $@ $(PREFIX)/bin/; \
	else \
		$(ECHO) $(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $@.o $(LIBS); \
		$(CXX) $(LDFLAGS) $(LDFLAGS_GUI_MINGW) -o $@ $@.o $(LIBS); \
		$(ECHO) $(CP) $@ $(PREFIX)/bin/; \
		$(CP) $@ $(PREFIX)/bin/; \
	fi;


%: %.o
	@if [ -f $@_icon.ico -a "$(OS)" = "MINGW32" ]; then \
		$(ECHO) "Windows & Icon detected ... embedding icon into binary."; \
		$(MAKE) $@_icon.o; \
		$(ECHO) $(CXX) $(LDFLAGS) -o $@ $< $@_icon.o $(LIBS); \
		$(CXX) $(LDFLAGS) -o $@ $< $@_icon.o $(LIBS); \
		$(ECHO) $(CP) $@ $(PREFIX)/bin/; \
		$(CP) $@ $(PREFIX)/bin/; \
	else \
		$(ECHO) $(CXX) $(LDFLAGS) -o $@ $< $(LIBS); \
		$(CXX) $(LDFLAGS) -o $@ $< $(LIBS); \
		$(ECHO) $(CP) $@ $(PREFIX)/bin/; \
		$(CP) $@ $(PREFIX)/bin/; \
	fi;



# For debugging the TOPP applications, please create the file 'Makefile_debug'
# in source/APPLICATIONS/ and add the line 'ADD_CXXFLAGS=-DDEBUG_TOPP' to it
-include ../Makefile_debug
-include Makefile_debug
