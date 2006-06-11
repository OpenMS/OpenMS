
# -*- Mode: make; tab-width: 2; -*-
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
# $Id: Makefile,v 1.25 2006/06/08 14:29:18 marc_sturm Exp $
# $Author: marc_sturm $
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include ../config_defs.mak

DIRECTORY=VISUAL

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

# Where the UIC-created headers reside
UIC_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)/UIC


THISLIB=$(LIBNAME)

# These depend on QT:
GUI_CPP_SOURCES=\
	BinnedRepWidget.C\
	IgnoringCombo.C\
	SpectrumCanvas.C\
	SpectrumWidget.C\
	SpectrumWindow.C\
	Spectrum1DCanvas.C\
	Spectrum1DWidget.C\
	Spectrum1DWindow.C\
	Spectrum2DCanvas.C\
	Spectrum2DWidget.C\
	Spectrum2DWindow.C\
	Spectrum3DCanvas.C\
	Spectrum3DOpenGLCanvas.C\
	Spectrum3DWidget.C\
	Spectrum3DWindow.C\
	SpectrumMDIWindow.C\
	EnhancedTabBar.C\
	AxisWidget.C\
	ColorSelector.C\
	MultiGradientSelector.C\
	MultiGradient.C\
	HistogramWidget.C\
	PreferencesManager.C\
	EnhancedWidgetStack.C\
	ListStack.C\
	PeakIcon.C\
	MappingInfo.C\
	FactoryProductView.C\
	SimMatrixWidget.C\
	ScaleWidget.C\
	ResultView.C\
	ClusterExperimentView.C\
	ClusterAnalysisWindow.C\
	ClusterRunWidget.C\
	LayerItem.C\
	LayerManager.C\
	AxisTickCalculator.C		

GUI_MOC_SOURCES=\
	BinnedRepWidget_moc.C\
	SpectrumCanvas_moc.C\
	SpectrumWindow_moc.C\
	SpectrumWidget_moc.C\
	Spectrum1DCanvas_moc.C\
	Spectrum1DWidget_moc.C\
	Spectrum1DWindow_moc.C\
	Spectrum2DCanvas_moc.C\
	Spectrum2DWidget_moc.C\
	Spectrum2DWindow_moc.C\
	Spectrum3DCanvas_moc.C\
	Spectrum3DOpenGLCanvas_moc.C\
	Spectrum3DWidget_moc.C\
	Spectrum3DWindow_moc.C\
	SpectrumMDIWindow_moc.C\
	EnhancedTabBar_moc.C\
	AxisWidget_moc.C\
	ColorSelector_moc.C\
	MultiGradientSelector_moc.C\
	HistogramWidget_moc.C\
	EnhancedWidgetStack_moc.C\
	ListStack_moc.C\
	FactoryProductView_moc.C\
	SimMatrixWidget_moc.C\
	ScaleWidget_moc.C\
	ResultView_moc.C\
	ClusterExperimentView_moc.C\
	ClusterAnalysisWindow_moc.C\
	ClusterRunWidget_moc.C\
	LayerItemTemplate_moc.C\
	LayerItem_moc.C\
	LayerManagerTemplate_moc.C\
	LayerManager_moc.C
	
UIC_SOURCES=\
	LayerItemTemplate.C\
	LayerManagerTemplate.C

#conditional compiles QT
ifeq ($(ENABLE_GUI),true)
	ADD_CPP_SOURCES=$(GUI_CPP_SOURCES)
	ADD_MOC_SOURCES=$(GUI_MOC_SOURCES)
endif

ADD_INCLUDES=$(GUI_INCLUDES)

SOURCES		= ${UIC_SOURCES} $(CPP_SOURCES) $(MOC_SOURCES) $(ADD_CPP_SOURCES) $(ADD_MOC_SOURCES) $(ANDI_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

# Take care of subdirectories of VISUAL

VISUAL_DIRS=DIALOGS

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

include ../config.mak

BinnedRepWidget_moc.C:  $(HEADER_DIR)/BinnedRepWidget.h
	$(MOC) $(HEADER_DIR)/BinnedRepWidget.h -o BinnedRepWidget_moc.C

ClusterAnalysisWindow_moc.C:  $(HEADER_DIR)/ClusterAnalysisWindow.h
	$(MOC) $(HEADER_DIR)/ClusterAnalysisWindow.h -o ClusterAnalysisWindow_moc.C

ClusterExperimentView_moc.C:  $(HEADER_DIR)/ClusterExperimentView.h
	$(MOC) $(HEADER_DIR)/ClusterExperimentView.h -o ClusterExperimentView_moc.C

ClusterRunWidget_moc.C:  $(HEADER_DIR)/ClusterRunWidget.h
	$(MOC) $(HEADER_DIR)/ClusterRunWidget.h -o ClusterRunWidget_moc.C

FactoryProductView_moc.C:  $(HEADER_DIR)/FactoryProductView.h
	$(MOC) $(HEADER_DIR)/FactoryProductView.h -o FactoryProductView_moc.C

ResultView_moc.C:  $(HEADER_DIR)/ResultView.h
	$(MOC) $(HEADER_DIR)/ResultView.h -o ResultView_moc.C

ScaleWidget_moc.C:  $(HEADER_DIR)/ScaleWidget.h
	$(MOC) $(HEADER_DIR)/ScaleWidget.h -o ScaleWidget_moc.C

SimMatrixWidget_moc.C:  $(HEADER_DIR)/SimMatrixWidget.h
	$(MOC) $(HEADER_DIR)/SimMatrixWidget.h -o SimMatrixWidget_moc.C

SpectrumCanvas_moc.C:  $(HEADER_DIR)/SpectrumCanvas.h
	$(MOC) $(HEADER_DIR)/SpectrumCanvas.h -o SpectrumCanvas_moc.C

SpectrumWindow_moc.C:  $(HEADER_DIR)/SpectrumWindow.h
	$(MOC) $(HEADER_DIR)/SpectrumWindow.h -o SpectrumWindow_moc.C

SpectrumWidget_moc.C:  $(HEADER_DIR)/SpectrumWidget.h
	$(MOC) $(HEADER_DIR)/SpectrumWidget.h -o SpectrumWidget_moc.C

Spectrum1DCanvas_moc.C:  $(HEADER_DIR)/Spectrum1DCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum1DCanvas.h -o Spectrum1DCanvas_moc.C

Spectrum1DWidget_moc.C:  $(HEADER_DIR)/Spectrum1DWidget.h
	$(MOC) $(HEADER_DIR)/Spectrum1DWidget.h -o Spectrum1DWidget_moc.C

Spectrum1DWindow_moc.C:  $(HEADER_DIR)/Spectrum1DWindow.h
	$(MOC) $(HEADER_DIR)/Spectrum1DWindow.h -o Spectrum1DWindow_moc.C

Spectrum2DCanvas_moc.C:  $(HEADER_DIR)/Spectrum2DCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum2DCanvas.h -o Spectrum2DCanvas_moc.C

Spectrum2DWidget_moc.C:  $(HEADER_DIR)/Spectrum2DWidget.h
	$(MOC) $(HEADER_DIR)/Spectrum2DWidget.h -o Spectrum2DWidget_moc.C

Spectrum2DWindow_moc.C:  $(HEADER_DIR)/Spectrum2DWindow.h
	$(MOC) $(HEADER_DIR)/Spectrum2DWindow.h -o Spectrum2DWindow_moc.C

Spectrum3DCanvas_moc.C:  $(HEADER_DIR)/Spectrum3DCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum3DCanvas.h -o Spectrum3DCanvas_moc.C

Spectrum3DOpenGLCanvas_moc.C:  $(HEADER_DIR)/Spectrum3DOpenGLCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum3DOpenGLCanvas.h -o Spectrum3DOpenGLCanvas_moc.C

Spectrum3DWidget_moc.C:  $(HEADER_DIR)/Spectrum3DWidget.h
	$(MOC) $(HEADER_DIR)/Spectrum3DWidget.h -o Spectrum3DWidget_moc.C

Spectrum3DWindow_moc.C:  $(HEADER_DIR)/Spectrum3DWindow.h
	$(MOC) $(HEADER_DIR)/Spectrum3DWindow.h -o Spectrum3DWindow_moc.C

SpectrumMDIWindow_moc.C:  $(HEADER_DIR)/SpectrumMDIWindow.h
	$(MOC) $(HEADER_DIR)/SpectrumMDIWindow.h -o SpectrumMDIWindow_moc.C

EnhancedTabBar_moc.C:  $(HEADER_DIR)/EnhancedTabBar.h
	$(MOC) $(HEADER_DIR)/EnhancedTabBar.h -o EnhancedTabBar_moc.C

AxisWidget_moc.C:  $(HEADER_DIR)/AxisWidget.h
	$(MOC) $(HEADER_DIR)/AxisWidget.h -o AxisWidget_moc.C

ColorSelector_moc.C:  $(HEADER_DIR)/ColorSelector.h
	$(MOC) $(HEADER_DIR)/ColorSelector.h -o ColorSelector_moc.C

MultiGradientSelector_moc.C:  $(HEADER_DIR)/MultiGradientSelector.h
	$(MOC) $(HEADER_DIR)/MultiGradientSelector.h -o MultiGradientSelector_moc.C

HistogramWidget_moc.C:  $(HEADER_DIR)/HistogramWidget.h
	$(MOC) $(HEADER_DIR)/HistogramWidget.h -o HistogramWidget_moc.C

EnhancedWidgetStack_moc.C:  $(HEADER_DIR)/EnhancedWidgetStack.h
	$(MOC) $(HEADER_DIR)/EnhancedWidgetStack.h -o EnhancedWidgetStack_moc.C	

ListStack_moc.C:  $(HEADER_DIR)/ListStack.h
	$(MOC) $(HEADER_DIR)/ListStack.h -o ListStack_moc.C	

LayerItem_moc.C:  $(HEADER_DIR)/LayerItem.h
	$(MOC) $(HEADER_DIR)/LayerItem.h -o LayerItem_moc.C

LayerManager_moc.C:  $(HEADER_DIR)/LayerManager.h
	$(MOC) $(HEADER_DIR)/LayerManager.h -o LayerManager_moc.C

