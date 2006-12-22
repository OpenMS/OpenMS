# -*- Mode: make; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
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

DIRECTORY=VISUAL/DIALOGS

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

# Where the UIC-created headers reside
UIC_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)/UIC

THISLIB=$(LIBNAME)

# These depend on QT:
GUI_CPP_SOURCES=\
	PreferencesDialog.C\
	PreferencesDialogPage.C\
	TOPPViewBasePDP.C\
	Spectrum1DWidgetPDP.C\
	Spectrum1DCanvasPDP.C\
	Spectrum2DWidgetPDP.C\
	Spectrum2DCanvasPDP.C\
	Spectrum3DWidgetPDP.C\
  Spectrum3DCanvasPDP.C\
	DBDialog.C\
	SaveImageDialog.C\
	DBSpectrumSelectorDialog.C\
	HistogramDialog.C\
	FactoryProductDialog.C\
	InspectDialog.C\
	Spectrum1DGoToDialog.C\
	Spectrum2DGoToDialog.C\
	OpenDialog.C\
	PeakPickingDialog.C\
	SmoothingDialog.C\
	BaselineFilteringDialog.C\
	SortDialog.C\
	FeaFiDialog.C


GUI_MOC_SOURCES=\
	PreferencesDialog_moc.C\
	PreferencesDialogPage_moc.C\
	Spectrum1DWidgetPDP_moc.C\
	Spectrum1DCanvasPDP_moc.C\
	Spectrum2DWidgetPDP_moc.C\
	Spectrum2DCanvasPDP_moc.C\
	Spectrum3DWidgetPDP_moc.C\
	Spectrum3DCanvasPDP_moc.C\
	DBDialog_moc.C\
	SaveImageDialog_moc.C\
	DBSpectrumSelectorDialog_moc.C\
	HistogramDialog_moc.C\
	FactoryProductDialog_moc.C\
	InspectDialog_moc.C\
	Spectrum1DGoToDialogTemplate_moc.C\
	Spectrum1DGoToDialog_moc.C\
	Spectrum2DGoToDialogTemplate_moc.C\
	Spectrum2DGoToDialog_moc.C\
	OpenDialogTemplate_moc.C\
	OpenDialog_moc.C\
	PeakPickingDialogTemplate_moc.C\
	PeakPickingDialog_moc.C\
	SmoothingDialogTemplate_moc.C\
	SmoothingDialog_moc.C\
	BaselineFilteringDialogTemplate_moc.C\
	BaselineFilteringDialog_moc.C\
	SortDialog_moc.C\
	FeaFiDialogTemplate_moc.C\
	FeaFiDialog_moc.C\
	TOPPViewBasePDP_moc.C

GUI_UIC_SOURCES=\
	Spectrum1DGoToDialogTemplate.C\
	Spectrum2DGoToDialogTemplate.C\
	OpenDialogTemplate.C\
	PeakPickingDialogTemplate.C\
	SmoothingDialogTemplate.C\
	BaselineFilteringDialogTemplate.C\
	FeaFiDialogTemplate.C

#conditional compiles QT
ifeq ($(ENABLE_GUI),true)
	ADD_CPP_SOURCES=\
		$(GUI_CPP_SOURCES)
	ADD_MOC_SOURCES=\
		$(GUI_MOC_SOURCES)
	UIC_SOURCES=\
		$(GUI_UIC_SOURCES)
endif

ADD_INCLUDES=$(GUI_INCLUDES)

SOURCES		= ${UIC_SOURCES} $(CPP_SOURCES) $(MOC_SOURCES) $(ADD_CPP_SOURCES) $(ADD_MOC_SOURCES) $(ANDI_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

-include .Dependencies

default: echodir $(OBJECTS)

depend:	$(UIC_SOURCES) $(MOC_SOURCES) $(ADD_MOC_SOURCES)
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX) -M  -fPIC -pedantic -Wno-long-long -Wno-deprecated  $(CPP_SOURCES) $(ADD_CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies

uic_depend: $(UIC_SOURCES)

uic_clean:
	@echo "Cleaning UIC files in $(UIC_DIR) and local direcory..."
	@$(RM) -f $(UIC_DIR)/*.h
	@$(RM) -f *Template.C

default_clean: uic_clean
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES) $(ADD_MOC_SOURCES)

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/$(THISLIB).objects ; done

include ../../config.mak

TOPPViewBasePDP_moc.C:  $(HEADER_DIR)/TOPPViewBasePDP.h
	$(MOC) $(HEADER_DIR)/TOPPViewBasePDP.h -o TOPPViewBasePDP_moc.C

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

DBDialog_moc.C:  $(HEADER_DIR)/DBDialog.h
	$(MOC) $(HEADER_DIR)/DBDialog.h -o DBDialog_moc.C

InspectDialog_moc.C:  $(HEADER_DIR)/InspectDialog.h
	$(MOC) $(HEADER_DIR)/InspectDialog.h -o InspectDialog_moc.C

SortDialog_moc.C:  $(HEADER_DIR)/SortDialog.h
	$(MOC) $(HEADER_DIR)/SortDialog.h -o SortDialog_moc.C

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

FeaFiDialog_moc.C:  $(HEADER_DIR)/FeaFiDialog.h
	$(MOC) $(HEADER_DIR)/FeaFiDialog.h -o FeaFiDialog_moc.C

Spectrum1DGoToDialog_moc.C:  $(HEADER_DIR)/Spectrum1DGoToDialog.h
	$(MOC) $(HEADER_DIR)/Spectrum1DGoToDialog.h -o Spectrum1DGoToDialog_moc.C

Spectrum2DGoToDialog_moc.C:  $(HEADER_DIR)/Spectrum2DGoToDialog.h
	$(MOC) $(HEADER_DIR)/Spectrum2DGoToDialog.h -o Spectrum2DGoToDialog_moc.C

OpenDialog_moc.C:  $(HEADER_DIR)/OpenDialog.h
	$(MOC) $(HEADER_DIR)/OpenDialog.h -o OpenDialog_moc.C

PeakPickingDialog_moc.C:  $(HEADER_DIR)/PeakPickingDialog.h
	$(MOC) $(HEADER_DIR)/PeakPickingDialog.h -o PeakPickingDialog_moc.C

SmoothingDialog_moc.C:  $(HEADER_DIR)/SmoothingDialog.h
	$(MOC) $(HEADER_DIR)/SmoothingDialog.h -o SmoothingDialog_moc.C

BaselineFilteringDialog_moc.C:  $(HEADER_DIR)/BaselineFilteringDialog.h
	$(MOC) $(HEADER_DIR)/BaselineFilteringDialog.h -o BaselineFilteringDialog_moc.C



