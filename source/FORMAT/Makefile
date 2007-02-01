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

THISLIB=$(LIBNAME)

-include .Dependencies

CPP_SOURCES=\
	Base64.C\
	DTAFile.C\
	DTA2DFile.C\
	PersistentObject.C\
	MascotInfile.C\
	MascotOutfile.C\
	InspectInfile.C\
	InspectOutfile.C\
	SequestInfile.C\
	SequestOutfile.C\
	TextFile.C\
	UniqueIdGenerator.C\
	FASTAFile.C\
	FileHandler.C\
	MzXMLFile.C\
	Param.C\
	MzDataFile.C\
	AnalysisXMLFile.C\
	DFeatureMapFile.C\
	DFeaturePairsFile.C\
	MascotXMLFile.C\
	DGridFile.C\
	SchemaFile.C\
	PeakFileOptions.C\
	PeakTypeEstimator.C\
	OMSSAXMLFile.C\
	OMSSACSVFile.C

ifeq ($(ANDIMS_SUPPORT),true)
	ifeq ($(NETCDF_SUPPORT),true)
		ANDI_SOURCES=\
			ANDIFile.C
	endif
endif

ifeq ($(LIBSVM_SUPPORT),true)
	LIBSVM_SOURCES=\
		LibSVMEncoder.C
endif

ifeq ($(CGAL_SUPPORT),true)
	CGAL_SOURCES=\
		ConsensusXMLFile.C
endif


# Take care of subdirectories of FORMAT

FORMAT_DIRS=HANDLERS DB

default: dirs $(OBJECTS)

dirs:	dummy
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET="default"

subdirs:	dummy
	@set -e; for i in $(SUBDIRS); do cd $$i; if $(MAKE) $(TARGET); then cd .. ; else cd .. ; exit 1 ; fi ; done

dummy:

depend:	dummy
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET=depend

default_clean:
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(EXECUTABLES) $(EXECUTABLE)
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET=default_clean

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/$(THISLIB).objects ; done
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET=collectlib


include ../config.mak

SOURCES		= $(CPP_SOURCES) $(ANDI_SOURCES) $(LIBSVM_SOURCES) $(CGAL_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

default:	echodir $(OBJECTS)

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."
