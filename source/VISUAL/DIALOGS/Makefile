# -*- Mode: make; tab-width: 2; -*-
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
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include ../../config_defs.mak

DIRECTORY=VISUAL/DIALOGS

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

# Where the UIC-created headers reside
UIC_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)/UIC


CPP_SOURCES=\
	PreferencesDialog.C\
	PreferencesDialogPage.C\
	Spectrum1DWidgetPDP.C\
	Spectrum1DCanvasPDP.C\
	Spectrum2DWidgetPDP.C\
	Spectrum2DCanvasPDP.C\
	Spectrum3DWidgetPDP.C\
  Spectrum3DCanvasPDP.C\
	SaveImageDialog.C\
	HistogramDialog.C\
	Spectrum1DGoToDialog.C\
	Spectrum2DGoToDialog.C\
	OpenDialog.C\
	DBSpectrumSelectorDialog.C\
	ToolsDialog.C\
	TOPPViewPrefDialog.C

MOC_SOURCES=\
	PreferencesDialog_moc.C\
	PreferencesDialogPage_moc.C\
	Spectrum1DWidgetPDP_moc.C\
	Spectrum1DCanvasPDP_moc.C\
	Spectrum2DWidgetPDP_moc.C\
	Spectrum2DCanvasPDP_moc.C\
	Spectrum3DWidgetPDP_moc.C\
	Spectrum3DCanvasPDP_moc.C\
	SaveImageDialog_moc.C\
	DBSpectrumSelectorDialog_moc.C\
	HistogramDialog_moc.C\
	Spectrum1DGoToDialog_moc.C\
	Spectrum2DGoToDialog_moc.C\
	OpenDialog_moc.C\
	ToolsDialog_moc.C\
	TOPPViewPrefDialog_moc.C

SOURCES		= $(CPP_SOURCES) $(MOC_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

-include .Dependencies

default: echodir $(OBJECTS)

depend: uic_depend $(MOC_SOURCES) 
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies

uic_depend: $(UIC_DIR)/Spectrum1DGoToDialogTemplate.h $(UIC_DIR)/Spectrum2DGoToDialogTemplate.h $(UIC_DIR)/OpenDialogTemplate.h $(UIC_DIR)/TOPPViewPrefDialogTemplate.h


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

AxisWidgetPDP_moc.C:  $(HEADER_DIR)/AxisWidgetPDP.h
	$(MOC) $(HEADER_DIR)/AxisWidgetPDP.h -o AxisWidgetPDP_moc.C

Spectrum1DCanvasPDP_moc.C:  $(HEADER_DIR)/Spectrum1DCanvasPDP.h
	$(MOC) $(HEADER_DIR)/Spectrum1DCanvasPDP.h -o Spectrum1DCanvasPDP_moc.C

Spectrum2DCanvasPDP_moc.C:  $(HEADER_DIR)/Spectrum2DCanvasPDP.h
	$(MOC) $(HEADER_DIR)/Spectrum2DCanvasPDP.h -o Spectrum2DCanvasPDP_moc.C

Spectrum3DCanvasPDP_moc.C:  $(HEADER_DIR)/Spectrum3DCanvasPDP.h
	$(MOC) $(HEADER_DIR)/Spectrum3DCanvasPDP.h -o Spectrum3DCanvasPDP_moc.C

Spectrum2DWidgetPDP_moc.C:  $(HEADER_DIR)/Spectrum2DWidgetPDP.h
	$(MOC) $(HEADER_DIR)/Spectrum2DWidgetPDP.h -o Spectrum2DWidgetPDP_moc.C

Spectrum1DWidgetPDP_moc.C:  $(HEADER_DIR)/Spectrum1DWidgetPDP.h
	$(MOC) $(HEADER_DIR)/Spectrum1DWidgetPDP.h -o Spectrum1DWidgetPDP_moc.C

Spectrum3DWidgetPDP_moc.C:  $(HEADER_DIR)/Spectrum3DWidgetPDP.h
	$(MOC) $(HEADER_DIR)/Spectrum3DWidgetPDP.h -o Spectrum3DWidgetPDP_moc.C

FactoryProductDialog_moc.C:  $(HEADER_DIR)/FactoryProductDialog.h
	$(MOC) $(HEADER_DIR)/FactoryProductDialog.h -o FactoryProductDialog_moc.C

SaveImageDialog_moc.C:  $(HEADER_DIR)/SaveImageDialog.h
	$(MOC) $(HEADER_DIR)/SaveImageDialog.h -o SaveImageDialog_moc.C

DBSpectrumSelectorDialog_moc.C:  $(HEADER_DIR)/DBSpectrumSelectorDialog.h
	$(MOC) $(HEADER_DIR)/DBSpectrumSelectorDialog.h -o DBSpectrumSelectorDialog_moc.C

HistogramDialog_moc.C:  $(HEADER_DIR)/HistogramDialog.h
	$(MOC) $(HEADER_DIR)/HistogramDialog.h -o HistogramDialog_moc.C

PreferencesDialog_moc.C:  $(HEADER_DIR)/PreferencesDialog.h
	$(MOC) $(HEADER_DIR)/PreferencesDialog.h -o PreferencesDialog_moc.C	

PreferencesDialogPage_moc.C:  $(HEADER_DIR)/PreferencesDialogPage.h
	$(MOC) $(HEADER_DIR)/PreferencesDialogPage.h -o PreferencesDialogPage_moc.C

ToolsDialog_moc.C:	$(HEADER_DIR)/ToolsDialog.h
	$(MOC) $(HEADER_DIR)/ToolsDialog.h -o ToolsDialog_moc.C

## UIC MOC SOURCES

TOPPViewPrefDialog_moc.C: $(HEADER_DIR)/TOPPViewPrefDialog.h $(UIC_DIR)/TOPPViewPrefDialogTemplate.h
	$(MOC) $(HEADER_DIR)/TOPPViewPrefDialog.h -o TOPPViewPrefDialog_moc.C

Spectrum1DGoToDialog_moc.C:  $(HEADER_DIR)/Spectrum1DGoToDialog.h
	$(MOC) $(HEADER_DIR)/Spectrum1DGoToDialog.h -o Spectrum1DGoToDialog_moc.C

Spectrum2DGoToDialog_moc.C:  $(HEADER_DIR)/Spectrum2DGoToDialog.h
	$(MOC) $(HEADER_DIR)/Spectrum2DGoToDialog.h -o Spectrum2DGoToDialog_moc.C

OpenDialog_moc.C:  $(HEADER_DIR)/OpenDialog.h
	$(MOC) $(HEADER_DIR)/OpenDialog.h -o OpenDialog_moc.C

## UIC HEADER

$(UIC_DIR)/TOPPViewPrefDialogTemplate.h: TOPPViewPrefDialog.ui
	$(UIC) -o $(UIC_DIR)/TOPPViewPrefDialogTemplate.h TOPPViewPrefDialog.ui

$(UIC_DIR)/Spectrum1DGoToDialogTemplate.h: Spectrum1DGoToDialog.ui
	$(UIC) -o $(UIC_DIR)/Spectrum1DGoToDialogTemplate.h Spectrum1DGoToDialog.ui

$(UIC_DIR)/Spectrum2DGoToDialogTemplate.h: Spectrum2DGoToDialog.ui
	$(UIC) -o $(UIC_DIR)/Spectrum2DGoToDialogTemplate.h Spectrum2DGoToDialog.ui

$(UIC_DIR)/OpenDialogTemplate.h: OpenDialog.ui
	$(UIC) -o $(UIC_DIR)/OpenDialogTemplate.h OpenDialog.ui
