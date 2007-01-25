# -*- Mode: make; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003 -- Oliver Kohlbacher, Ulf Leser, Knut Reinert
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
# $Id: Makefile,v 1.0 2006/03/15 15:36:03 stefan_heess Exp $
# $Author: stefan_heess $
# $Maintainer: Stefan Heess $
# --------------------------------------------------------------------------

include ../../config_defs.mak

DIRECTORY=VISUAL/VISUALIZER

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

# Where the UIC-created headers reside
UIC_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)/UIC


THISLIB=$(LIBNAME)

# These depend on QT:
GUI_CPP_SOURCES=\
	BaseVisualizer.C\
	SampleVisualizer.C\
	DigestionVisualizer.C\
	ModificationVisualizer.C\
	TaggingVisualizer.C\
	HPLCVisualizer.C\
	MetaInfoVisualizer.C\
	GradientVisualizer.C\
	SoftwareVisualizer.C\
	SourceFileVisualizer.C\
	ContactPersonVisualizer.C\
	InstrumentVisualizer.C\
	IonSourceVisualizer.C\
	IonDetectorVisualizer.C\
	MassAnalyzerVisualizer.C\
	ProcessingMethodVisualizer.C\
	ProteinIdentificationVisualizer.C\
	ProteinHitVisualizer.C\
	PeptideHitVisualizer.C\
	ExperimentalSettingsVisualizer.C\
	AcquisitionVisualizer.C\
	AcquisitionInfoVisualizer.C\
	MetaInfoDescriptionVisualizer.C\
	PrecursorVisualizer.C\
	InstrumentSettingsVisualizer.C\
	IdentificationVisualizer.C\
	SpectrumSettingsVisualizer.C\

GUI_MOC_SOURCES=\
	BaseVisualizer_moc.C\
	SampleVisualizer_moc.C\
	DigestionVisualizer_moc.C\
	ModificationVisualizer_moc.C\
	TaggingVisualizer_moc.C\
	HPLCVisualizer_moc.C\
	MetaInfoVisualizer_moc.C\
	GradientVisualizer_moc.C\
	SoftwareVisualizer_moc.C\
	SourceFileVisualizer_moc.C\
	ContactPersonVisualizer_moc.C\
	InstrumentVisualizer_moc.C\
	IonSourceVisualizer_moc.C\
	IonDetectorVisualizer_moc.C\
	MassAnalyzerVisualizer_moc.C\
  ProcessingMethodVisualizer_moc.C\
	ProteinIdentificationVisualizer_moc.C\
	ProteinHitVisualizer_moc.C\
	PeptideHitVisualizer_moc.C\
	ExperimentalSettingsVisualizer_moc.C\
	AcquisitionVisualizer_moc.C\
	AcquisitionInfoVisualizer_moc.C\
	MetaInfoDescriptionVisualizer_moc.C\
	PrecursorVisualizer_moc.C\
	InstrumentSettingsVisualizer_moc.C\
	IdentificationVisualizer_moc.C\
	SpectrumSettingsVisualizer_moc.C
	
UIC_SOURCES=\
	

#conditional compiles QT
ifeq ($(ENABLE_GUI),true)
	ADD_CPP_SOURCES=$(GUI_CPP_SOURCES)
	ADD_MOC_SOURCES=$(GUI_MOC_SOURCES)
endif

ADD_INCLUDES=$(GUI_INCLUDES)

SOURCES		= ${UIC_SOURCES} $(CPP_SOURCES) $(MOC_SOURCES) $(ADD_CPP_SOURCES) $(ADD_MOC_SOURCES) $(ANDI_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)


# Take care of subdirectories of VISUAL

VISUAL_DIRS=\
        
default:	dirs $(OBJECTS)

-include .Dependencies

dirs:	dummy
	@$(MAKE) subdirs SUBDIRS="$(VISUAL_DIRS)" TARGET="default"

subdirs:	dummy
	@set -e; for i in $(SUBDIRS); do cd $$i; if $(MAKE) $(TARGET); then cd .. ; else cd .. ; exit 1 ; fi ; done

dummy:

depend:	uic_depend $(MOC_SOURCES) $(ADD_MOC_SOURCES)
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies
	@$(MAKE) subdirs SUBDIRS="$(VISUAL_DIRS)" TARGET=depend

uic_depend: $(UIC_SOURCES)
	@$(MAKE) subdirs SUBDIRS="$(VISUAL_DIRS)" TARGET=uic_depend

uic_clean:
	@echo "Cleaning UIC files in $(UIC_DIR) and local direcory..."
	@$(RM) -f $(UIC_DIR)/*.h
	@$(RM) -f *Template.C

default_clean: uic_clean
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES) $(ADD_MOC_SOURCES)
	@$(MAKE) subdirs SUBDIRS="$(VISUAL_DIRS)" TARGET=default_clean

echodir:
	@echo "entering $(DIRECTORY)..."

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/$(THISLIB).objects ; done
	@$(MAKE) subdirs SUBDIRS="$(VISUAL_DIRS)" TARGET=collectlib

include ../../config.mak
	
BaseVisualizer_moc.C:  $(HEADER_DIR)/BaseVisualizer.h
	$(MOC) $(HEADER_DIR)/BaseVisualizer.h -o BaseVisualizer_moc.C	

SampleVisualizer_moc.C:  $(HEADER_DIR)/SampleVisualizer.h
	$(MOC) $(HEADER_DIR)/SampleVisualizer.h -o SampleVisualizer_moc.C
	
DigestionVisualizer_moc.C:  $(HEADER_DIR)/DigestionVisualizer.h
	$(MOC) $(HEADER_DIR)/DigestionVisualizer.h -o DigestionVisualizer_moc.C
	
ModificationVisualizer_moc.C:  $(HEADER_DIR)/ModificationVisualizer.h
	$(MOC) $(HEADER_DIR)/ModificationVisualizer.h -o ModificationVisualizer_moc.C
	
TaggingVisualizer_moc.C:  $(HEADER_DIR)/TaggingVisualizer.h
	$(MOC) $(HEADER_DIR)/TaggingVisualizer.h -o TaggingVisualizer_moc.C

HPLCVisualizer_moc.C:  $(HEADER_DIR)/HPLCVisualizer.h
	$(MOC) $(HEADER_DIR)/HPLCVisualizer.h -o HPLCVisualizer_moc.C

MetaInfoVisualizer_moc.C:  $(HEADER_DIR)/MetaInfoVisualizer.h
	$(MOC) $(HEADER_DIR)/MetaInfoVisualizer.h -o MetaInfoVisualizer_moc.C	

GradientVisualizer_moc.C:  $(HEADER_DIR)/GradientVisualizer.h
	$(MOC) $(HEADER_DIR)/GradientVisualizer.h -o GradientVisualizer_moc.C	

SoftwareVisualizer_moc.C:  $(HEADER_DIR)/SoftwareVisualizer.h
	$(MOC) $(HEADER_DIR)/SoftwareVisualizer.h -o SoftwareVisualizer_moc.C	

SourceFileVisualizer_moc.C:  $(HEADER_DIR)/SourceFileVisualizer.h
	$(MOC) $(HEADER_DIR)/SourceFileVisualizer.h -o SourceFileVisualizer_moc.C	

ContactPersonVisualizer_moc.C:  $(HEADER_DIR)/ContactPersonVisualizer.h
	$(MOC) $(HEADER_DIR)/ContactPersonVisualizer.h -o ContactPersonVisualizer_moc.C	

InstrumentVisualizer_moc.C:  $(HEADER_DIR)/InstrumentVisualizer.h
	$(MOC) $(HEADER_DIR)/InstrumentVisualizer.h -o InstrumentVisualizer_moc.C	

IonSourceVisualizer_moc.C:  $(HEADER_DIR)/IonSourceVisualizer.h
	$(MOC) $(HEADER_DIR)/IonSourceVisualizer.h -o IonSourceVisualizer_moc.C	

IonDetectorVisualizer_moc.C:  $(HEADER_DIR)/IonDetectorVisualizer.h
	$(MOC) $(HEADER_DIR)/IonDetectorVisualizer.h -o IonDetectorVisualizer_moc.C		
	
MassAnalyzerVisualizer_moc.C:  $(HEADER_DIR)/MassAnalyzerVisualizer.h
	$(MOC) $(HEADER_DIR)/MassAnalyzerVisualizer.h -o MassAnalyzerVisualizer_moc.C		
	
ProcessingMethodVisualizer_moc.C:  $(HEADER_DIR)/ProcessingMethodVisualizer.h
	$(MOC) $(HEADER_DIR)/ProcessingMethodVisualizer.h -o ProcessingMethodVisualizer_moc.C		

ProteinIdentificationVisualizer_moc.C:  $(HEADER_DIR)/ProteinIdentificationVisualizer.h
	$(MOC) $(HEADER_DIR)/ProteinIdentificationVisualizer.h -o ProteinIdentificationVisualizer_moc.C		

ProteinHitVisualizer_moc.C:  $(HEADER_DIR)/ProteinHitVisualizer.h
	$(MOC) $(HEADER_DIR)/ProteinHitVisualizer.h -o ProteinHitVisualizer_moc.C		

PeptideHitVisualizer_moc.C:  $(HEADER_DIR)/PeptideHitVisualizer.h
	$(MOC) $(HEADER_DIR)/PeptideHitVisualizer.h -o PeptideHitVisualizer_moc.C			
	
ExperimentalSettingsVisualizer_moc.C:  $(HEADER_DIR)/ExperimentalSettingsVisualizer.h
	$(MOC) $(HEADER_DIR)/ExperimentalSettingsVisualizer.h -o ExperimentalSettingsVisualizer_moc.C	

AcquisitionVisualizer_moc.C:  $(HEADER_DIR)/AcquisitionVisualizer.h
	$(MOC) $(HEADER_DIR)/AcquisitionVisualizer.h -o AcquisitionVisualizer_moc.C		
	
AcquisitionInfoVisualizer_moc.C:  $(HEADER_DIR)/AcquisitionInfoVisualizer.h
	$(MOC) $(HEADER_DIR)/AcquisitionInfoVisualizer.h -o AcquisitionInfoVisualizer_moc.C		

MetaInfoDescriptionVisualizer_moc.C:  $(HEADER_DIR)/MetaInfoDescriptionVisualizer.h
	$(MOC) $(HEADER_DIR)/MetaInfoDescriptionVisualizer.h -o MetaInfoDescriptionVisualizer_moc.C		

PrecursorVisualizer_moc.C:  $(HEADER_DIR)/PrecursorVisualizer.h
	$(MOC) $(HEADER_DIR)/PrecursorVisualizer.h -o PrecursorVisualizer_moc.C		

InstrumentSettingsVisualizer_moc.C:  $(HEADER_DIR)/InstrumentSettingsVisualizer.h
	$(MOC) $(HEADER_DIR)/InstrumentSettingsVisualizer.h -o InstrumentSettingsVisualizer_moc.C		

IdentificationVisualizer_moc.C:  $(HEADER_DIR)/IdentificationVisualizer.h
	$(MOC) $(HEADER_DIR)/IdentificationVisualizer.h -o IdentificationVisualizer_moc.C		
	
SpectrumSettingsVisualizer_moc.C:  $(HEADER_DIR)/SpectrumSettingsVisualizer.h
	$(MOC) $(HEADER_DIR)/SpectrumSettingsVisualizer.h -o SpectrumSettingsVisualizer_moc.C	