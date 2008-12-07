# -*- mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
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

DIRECTORY=VISUAL/DIALOGS

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

# Where the UIC-created headers reside
UIC_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)/UIC


CPP_SOURCES=\
	SaveImageDialog.C\
	HistogramDialog.C\
	Spectrum1DGoToDialog.C\
	Spectrum2DGoToDialog.C\
	TOPPViewOpenDialog.C\
	DBOpenDialog.C\
	ToolsDialog.C\
	TOPPViewPrefDialog.C\
	Spectrum1DPrefDialog.C\
	Spectrum2DPrefDialog.C\
	Spectrum3DPrefDialog.C\
	DataFilterDialog.C\
	LayerStatisticsDialog.C\
	TheoreticalSpectrumGenerationDialog.C\
	SpectrumAlignmentDialog.C\
	FeatureEditDialog.C

MOC_SOURCES=\
	SaveImageDialog_moc.C\
	DBOpenDialog_moc.C\
	HistogramDialog_moc.C\
	Spectrum1DGoToDialog_moc.C\
	Spectrum2DGoToDialog_moc.C\
	TOPPViewOpenDialog_moc.C\
	ToolsDialog_moc.C\
	TOPPViewPrefDialog_moc.C\
	Spectrum1DPrefDialog_moc.C\
	Spectrum2DPrefDialog_moc.C\
	Spectrum3DPrefDialog_moc.C\
	DataFilterDialog_moc.C\
	LayerStatisticsDialog_moc.C\
	TheoreticalSpectrumGenerationDialog_moc.C\
	SpectrumAlignmentDialog_moc.C\
	FeatureEditDialog_moc.C

SOURCES		= $(CPP_SOURCES) $(MOC_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

-include .Dependencies

default: echodir uic_depend $(OBJECTS)

depend: uic_depend $(MOC_SOURCES) 
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX_MAKEDEPEND) $(MAKEDEP_OPTS) $(CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies

uic_depend: $(UIC_DIR)/Spectrum1DGoToDialogTemplate.h $(UIC_DIR)/Spectrum2DGoToDialogTemplate.h $(UIC_DIR)/TOPPViewOpenDialogTemplate.h $(UIC_DIR)/TOPPViewPrefDialogTemplate.h $(UIC_DIR)/Spectrum1DPrefDialogTemplate.h $(UIC_DIR)/Spectrum2DPrefDialogTemplate.h $(UIC_DIR)/Spectrum3DPrefDialogTemplate.h $(UIC_DIR)/DataFilterDialogTemplate.h $(UIC_DIR)/LayerStatisticsDialogTemplate.h $(UIC_DIR)/TheoreticalSpectrumGenerationDialogTemplate.h $(UIC_DIR)/SpectrumAlignmentDialogTemplate.h $(UIC_DIR)/FeatureEditDialogTemplate.h

uic_clean:
	@echo "Cleaning UIC files in $(UIC_DIR) and local direcory..."
	@$(RM) -f $(UIC_DIR)/*.h

default_clean: uic_clean
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES)

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/libOpenMS.objects ; done

include ../../config.mak

## MOC SOURCES

SaveImageDialog_moc.C:  $(HEADER_DIR)/SaveImageDialog.h
	$(MOC) $(HEADER_DIR)/SaveImageDialog.h -o SaveImageDialog_moc.C

DBOpenDialog_moc.C:  $(HEADER_DIR)/DBOpenDialog.h
	$(MOC) $(HEADER_DIR)/DBOpenDialog.h -o DBOpenDialog_moc.C

HistogramDialog_moc.C:  $(HEADER_DIR)/HistogramDialog.h
	$(MOC) $(HEADER_DIR)/HistogramDialog.h -o HistogramDialog_moc.C

ToolsDialog_moc.C:	$(HEADER_DIR)/ToolsDialog.h
	$(MOC) $(HEADER_DIR)/ToolsDialog.h -o ToolsDialog_moc.C

## UIC MOC SOURCES

TOPPViewPrefDialog_moc.C: $(HEADER_DIR)/TOPPViewPrefDialog.h $(UIC_DIR)/TOPPViewPrefDialogTemplate.h
	$(MOC) $(HEADER_DIR)/TOPPViewPrefDialog.h -o TOPPViewPrefDialog_moc.C

Spectrum1DPrefDialog_moc.C: $(HEADER_DIR)/Spectrum1DPrefDialog.h $(UIC_DIR)/Spectrum1DPrefDialogTemplate.h
	$(MOC) $(HEADER_DIR)/Spectrum1DPrefDialog.h -o Spectrum1DPrefDialog_moc.C

Spectrum2DPrefDialog_moc.C: $(HEADER_DIR)/Spectrum2DPrefDialog.h $(UIC_DIR)/Spectrum2DPrefDialogTemplate.h
	$(MOC) $(HEADER_DIR)/Spectrum2DPrefDialog.h -o Spectrum2DPrefDialog_moc.C

Spectrum3DPrefDialog_moc.C: $(HEADER_DIR)/Spectrum3DPrefDialog.h $(UIC_DIR)/Spectrum3DPrefDialogTemplate.h
	$(MOC) $(HEADER_DIR)/Spectrum3DPrefDialog.h -o Spectrum3DPrefDialog_moc.C

Spectrum1DGoToDialog_moc.C:  $(HEADER_DIR)/Spectrum1DGoToDialog.h
	$(MOC) $(HEADER_DIR)/Spectrum1DGoToDialog.h -o Spectrum1DGoToDialog_moc.C

Spectrum2DGoToDialog_moc.C:  $(HEADER_DIR)/Spectrum2DGoToDialog.h
	$(MOC) $(HEADER_DIR)/Spectrum2DGoToDialog.h -o Spectrum2DGoToDialog_moc.C

TOPPViewOpenDialog_moc.C:  $(HEADER_DIR)/TOPPViewOpenDialog.h
	$(MOC) $(HEADER_DIR)/TOPPViewOpenDialog.h -o TOPPViewOpenDialog_moc.C

DataFilterDialog_moc.C:  $(HEADER_DIR)/DataFilterDialog.h
	$(MOC) $(HEADER_DIR)/DataFilterDialog.h -o DataFilterDialog_moc.C

LayerStatisticsDialog_moc.C:  $(HEADER_DIR)/LayerStatisticsDialog.h
	$(MOC) $(HEADER_DIR)/LayerStatisticsDialog.h -o LayerStatisticsDialog_moc.C

TheoreticalSpectrumGenerationDialog_moc.C:	$(HEADER_DIR)/TheoreticalSpectrumGenerationDialog.h
	$(MOC) $(HEADER_DIR)/TheoreticalSpectrumGenerationDialog.h -o TheoreticalSpectrumGenerationDialog_moc.C
	
SpectrumAlignmentDialog_moc.C:	$(HEADER_DIR)/SpectrumAlignmentDialog.h
	$(MOC) $(HEADER_DIR)/SpectrumAlignmentDialog.h -o SpectrumAlignmentDialog_moc.C

FeatureEditDialog_moc.C:  $(HEADER_DIR)/FeatureEditDialog.h
	$(MOC) $(HEADER_DIR)/FeatureEditDialog.h -o FeatureEditDialog_moc.C

## UIC HEADER

$(UIC_DIR)/TOPPViewPrefDialogTemplate.h: TOPPViewPrefDialog.ui
	$(UIC) -o $(UIC_DIR)/TOPPViewPrefDialogTemplate.h TOPPViewPrefDialog.ui

$(UIC_DIR)/Spectrum1DPrefDialogTemplate.h: Spectrum1DPrefDialog.ui
	$(UIC) -o $(UIC_DIR)/Spectrum1DPrefDialogTemplate.h Spectrum1DPrefDialog.ui

$(UIC_DIR)/Spectrum2DPrefDialogTemplate.h: Spectrum2DPrefDialog.ui
	$(UIC) -o $(UIC_DIR)/Spectrum2DPrefDialogTemplate.h Spectrum2DPrefDialog.ui

$(UIC_DIR)/Spectrum3DPrefDialogTemplate.h: Spectrum3DPrefDialog.ui
	$(UIC) -o $(UIC_DIR)/Spectrum3DPrefDialogTemplate.h Spectrum3DPrefDialog.ui

$(UIC_DIR)/Spectrum1DGoToDialogTemplate.h: Spectrum1DGoToDialog.ui
	$(UIC) -o $(UIC_DIR)/Spectrum1DGoToDialogTemplate.h Spectrum1DGoToDialog.ui

$(UIC_DIR)/Spectrum2DGoToDialogTemplate.h: Spectrum2DGoToDialog.ui
	$(UIC) -o $(UIC_DIR)/Spectrum2DGoToDialogTemplate.h Spectrum2DGoToDialog.ui

$(UIC_DIR)/TOPPViewOpenDialogTemplate.h: TOPPViewOpenDialog.ui
	$(UIC) -o $(UIC_DIR)/TOPPViewOpenDialogTemplate.h TOPPViewOpenDialog.ui

$(UIC_DIR)/DataFilterDialogTemplate.h: DataFilterDialog.ui
	$(UIC) -o $(UIC_DIR)/DataFilterDialogTemplate.h DataFilterDialog.ui

$(UIC_DIR)/LayerStatisticsDialogTemplate.h: LayerStatisticsDialog.ui
	$(UIC) -o $(UIC_DIR)/LayerStatisticsDialogTemplate.h LayerStatisticsDialog.ui

$(UIC_DIR)/TheoreticalSpectrumGenerationDialogTemplate.h: TheoreticalSpectrumGenerationDialog.ui
	$(UIC) -o $(UIC_DIR)/TheoreticalSpectrumGenerationDialogTemplate.h TheoreticalSpectrumGenerationDialog.ui
	
$(UIC_DIR)/SpectrumAlignmentDialogTemplate.h: SpectrumAlignmentDialog.ui
	$(UIC) -o $(UIC_DIR)/SpectrumAlignmentDialogTemplate.h SpectrumAlignmentDialog.ui

$(UIC_DIR)/FeatureEditDialogTemplate.h: FeatureEditDialog.ui
	$(UIC) -o $(UIC_DIR)/FeatureEditDialogTemplate.h FeatureEditDialog.ui

