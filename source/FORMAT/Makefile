# -*- Mode: makefile; tab-width: 2; -*-
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
# $Maintainer:  Early Pearly$
# --------------------------------------------------------------------------

all:	default

include ../config_defs.mak

DIRECTORY=FORMAT


-include .Dependencies

CPP_SOURCES=\
	PeakFileOptions.C\
	PeakTypeEstimator.C\
	UniqueIdGenerator.C\
	Base64.C\
	DTAFile.C\
	DTA2DFile.C\
	PersistentObject.C\
	TextFile.C\
	XMLFile.C\
	Schema2CV.C\
	ControlledVocabulary.C\
	MzXMLFile.C\
	MascotInfile.C\
	MascotOutfile.C\
	InspectInfile.C\
	InspectOutfile.C\
	SequestInfile.C\
	SequestOutfile.C\
	PepNovoInfile.C\
	PepNovoOutfile.C\
	FASTAFile.C\
	MzDataFile.C\
	IdXMLFile.C\
	FeatureXMLFile.C\
	FeaturePairsXMLFile.C\
	MascotXMLFile.C\
	GridFile.C\
	OMSSAXMLFile.C\
	OMSSACSVFile.C\
	PTMXMLFile.C\
	ConsensusXMLFile.C\
	MSPFile.C\
	XTandemXMLFile.C\
	XMLValidator.C\
	XTandemInfile.C\
	FastaIterator.C\
	FastaIteratorIntern.C\
	FileHandler.C

ifeq ($(ANDIMS_SUPPORT),true)
	ANDI_SOURCES=\
		ANDIFile.C
endif

ifeq ($(LIBSVM_SUPPORT),true)
	LIBSVM_SOURCES=\
		LibSVMEncoder.C
endif


# Redefine targets to take care of subdirectories of FORMAT

dummy:

depend:	dummy
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies
	@cd DB; if $(MAKE) depend; then cd .. ; else cd .. ; exit 1; fi;
	@cd HANDLERS; if $(MAKE) depend; then cd .. ; else cd .. ; exit 1; fi;

default_clean: dummy
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(EXECUTABLES) $(EXECUTABLE)
	@cd DB; if $(MAKE) default_clean; then cd .. ; else cd .. ; exit 1; fi;
	@cd HANDLERS; if $(MAKE) default_clean; then cd .. ; else cd .. ; exit 1; fi;

collectlib:	dummy
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/libOpenMS.objects ; done
	@cd DB; if $(MAKE) collectlib; then cd .. ; else cd .. ; exit 1; fi;
	@cd HANDLERS; if $(MAKE) collectlib; then cd .. ; else cd .. ; exit 1; fi;


include ../config.mak

SOURCES		= $(CPP_SOURCES) $(ANDI_SOURCES) $(LIBSVM_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

default: $(OBJECTS) dummy
	@echo "Entering directory $(DIRECTORY) ..."
	@cd DB; if $(MAKE) default; then cd .. ; else cd .. ; exit 1; fi;
	@cd HANDLERS; if $(MAKE) default; then cd .. ; else cd .. ; exit 1; fi;


