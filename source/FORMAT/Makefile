# -*- Mode: makefile; tab-width: 2; -*-
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
	InspectFile.C\
	InspectInfile.C\
	InspectOutfile.C\
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
	ConsensusXMLFile.C\
	SequestInfile.C\
	SequestOutfile.C

ifeq ($(ANDIMS_SUPPORT),true)
	ifeq ($(NETCDF_SUPPORT),true)
		ANDI_SOURCES=\
			ANDIFile.C
	endif
endif

ifeq ($(ENABLE_DB),true)
	DB_CPP_SOURCES=\
		DBConnection.C\
		DataSetInfo.C\
		DBAdapter.C
endif

ifeq ($(LIBSVM_SUPPORT),true)
	LIBSVM_SOURCES=\
		LibSVMEncoder.C
endif

# Take care of subdirectories of FORMAT

FORMAT_DIRS=HANDLERS

#conditional compiles QT
ifeq ($(ENABEL_DB),true)
	ADD_CPP_SOURCES= $(GUI_CPP_SOURCES) ${DB_CPP_SOURCES}
else
	ADD_CPP_SOURCES=  ${DB_CPP_SOURCES}
endif

default:	dirs $(OBJECTS)

dirs:	dummy
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET="default"

subdirs:	dummy
	@set -e; for i in $(SUBDIRS); do cd $$i; if $(MAKE) $(TARGET); then cd .. ; else cd .. ; exit 1 ; fi ; done

dummy:

depend:	dummy
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET=depend

default_clean:
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES) $(ADD_MOC_SOURCES) $(EXECUTABLES) $(EXECUTABLE)
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET=default_clean

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/$(THISLIB).objects ; done
	@$(MAKE) subdirs SUBDIRS="$(FORMAT_DIRS)" TARGET=collectlib


include ../config.mak

OPENMS_INCLUDES=$(GUI_INCLUDES)

SOURCES		= $(CPP_SOURCES) $(MOC_SOURCES) $(ADD_CPP_SOURCES) $(ADD_MOC_SOURCES) $(ANDI_SOURCES) $(LIBSVM_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

default:	echodir $(OBJECTS)

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."
