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

include ../config_defs.mak

DIRECTORY=VISUAL

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

# Where the UIC-created headers reside
UIC_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)/UIC



CPP_SOURCES=\
	EnhancedTabBar.C\
	AxisWidget.C\
	ColorSelector.C\
	MultiGradient.C\
	MultiGradientSelector.C\
	HistogramWidget.C\
	PeakIcon.C\
	AxisTickCalculator.C\
	MetaDataBrowser.C\
	LayerData.C\
	SpectrumCanvas.C\
	SpectrumWidget.C\
	Spectrum1DCanvas.C\
	Spectrum1DWidget.C\
	Spectrum2DCanvas.C\
	Spectrum2DWidget.C\
	Spectrum3DCanvas.C\
	Spectrum3DOpenGLCanvas.C\
	Spectrum3DWidget.C\
	ParamEditor.C\
	Annotations1DContainer.C\
	Annotation1DItem.C\
	Annotation1DDistanceItem.C\
	Annotation1DTextItem.C\
	Annotation1DPeakItem.C\


MOC_SOURCES=\
	EnhancedTabBar_moc.C\
	AxisWidget_moc.C\
	ColorSelector_moc.C\
	MultiGradientSelector_moc.C\
	HistogramWidget_moc.C\
	MetaDataBrowser_moc.C\
	SpectrumCanvas_moc.C\
	SpectrumWidget_moc.C\
	Spectrum1DCanvas_moc.C\
	Spectrum1DWidget_moc.C\
	Spectrum2DCanvas_moc.C\
	Spectrum2DWidget_moc.C\
	Spectrum3DCanvas_moc.C\
	Spectrum3DOpenGLCanvas_moc.C\
	Spectrum3DWidget_moc.C\
	ParamEditor_moc.C\

SOURCES		= $(CPP_SOURCES) $(MOC_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

# Redefine targets to take care of subdirectories of VISUAL

default: echodir uic_depend $(OBJECTS)
	@cd DIALOGS; if $(MAKE) default; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd VISUALIZER; if $(MAKE) default; then cd .. ; else cd .. ; exit 1 ; fi ;

-include .Dependencies

dummy:

depend:	uic_depend $(MOC_SOURCES)
	@echo "creating dependencies for $(DIRECTORY)..."
	@echo "" > .Dependencies
	@$(CXX_MAKEDEPEND) $(MAKEDEP_OPTS) $(MOC_SOURCES) $(CPP_SOURCES) $(STD_CPP_INCLUDES) $(ADD_INCLUDES) $(OPENMS_INCLUDES)  >.Dependencies
	@cd DIALOGS; if $(MAKE) depend; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd VISUALIZER; if $(MAKE) depend; then cd .. ; else cd .. ; exit 1 ; fi ;

uic_depend: $(UIC_DIR)/ParamEditorTemplate.h
	@cd DIALOGS; if $(MAKE) uic_depend; then cd .. ; else cd .. ; exit 1 ; fi ;

uic_clean:
	@echo "Cleaning UIC files in $(UIC_DIR) and local direcory..."
	@$(RM) -f $(UIC_DIR)/*.h
	@$(RM) -f *.dh
	@$(RM) -f *Template.C

default_clean: uic_clean
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES)
	@cd DIALOGS; if $(MAKE) default_clean; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd VISUALIZER; if $(MAKE) default_clean; then cd .. ; else cd .. ; exit 1 ; fi ;

echodir:
	@echo " "
	@echo "Entering directory $(DIRECTORY) ..."

collectlib:	echodir
	@echo "collecting objects from $(DIRECTORY)..."
	@for i in $(OBJECTS); do echo $(DIRECTORY)/$$i >>$(OPENMS_PATH)/source/libOpenMS.objects ; done
	@cd DIALOGS; if $(MAKE) collectlib; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd VISUALIZER; if $(MAKE) collectlib; then cd .. ; else cd .. ; exit 1 ; fi ;

include ../config.mak


SpectrumCanvas_moc.C:  $(HEADER_DIR)/SpectrumCanvas.h
	$(MOC) $(HEADER_DIR)/SpectrumCanvas.h -o SpectrumCanvas_moc.C

SpectrumWidget_moc.C:  $(HEADER_DIR)/SpectrumWidget.h
	$(MOC) $(HEADER_DIR)/SpectrumWidget.h -o SpectrumWidget_moc.C

Spectrum1DCanvas_moc.C:  $(HEADER_DIR)/Spectrum1DCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum1DCanvas.h -o Spectrum1DCanvas_moc.C

Spectrum1DWidget_moc.C:  $(HEADER_DIR)/Spectrum1DWidget.h
	$(MOC) $(HEADER_DIR)/Spectrum1DWidget.h -o Spectrum1DWidget_moc.C

Spectrum2DCanvas_moc.C:  $(HEADER_DIR)/Spectrum2DCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum2DCanvas.h -o Spectrum2DCanvas_moc.C

Spectrum2DWidget_moc.C:  $(HEADER_DIR)/Spectrum2DWidget.h
	$(MOC) $(HEADER_DIR)/Spectrum2DWidget.h -o Spectrum2DWidget_moc.C

Spectrum3DCanvas_moc.C:  $(HEADER_DIR)/Spectrum3DCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum3DCanvas.h -o Spectrum3DCanvas_moc.C

Spectrum3DOpenGLCanvas_moc.C:  $(HEADER_DIR)/Spectrum3DOpenGLCanvas.h
	$(MOC) $(HEADER_DIR)/Spectrum3DOpenGLCanvas.h -o Spectrum3DOpenGLCanvas_moc.C

Spectrum3DWidget_moc.C:  $(HEADER_DIR)/Spectrum3DWidget.h
	$(MOC) $(HEADER_DIR)/Spectrum3DWidget.h -o Spectrum3DWidget_moc.C

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

ListStack_moc.C:  $(HEADER_DIR)/ListStack.h
	$(MOC) $(HEADER_DIR)/ListStack.h -o ListStack_moc.C

MetaDataBrowser_moc.C:  $(HEADER_DIR)/MetaDataBrowser.h
	$(MOC) $(HEADER_DIR)/MetaDataBrowser.h -o MetaDataBrowser_moc.C

ParamEditor_moc.C:  $(HEADER_DIR)/ParamEditor.h
	$(MOC) $(HEADER_DIR)/ParamEditor.h -o ParamEditor_moc.C
	

# UIC headers


$(UIC_DIR)/ParamEditorTemplate.h: ParamEditor.ui
	$(UIC) -o $(UIC_DIR)/ParamEditorTemplate.h ParamEditor.ui



# For debugging the visual classes, create the file 'Makefile_debug' and 
# set 'ADD_CXXFLAGS' one or more of the following flags
#  -DDEBUG_TOPPVIEW  
#  -DTIMING_TOPPVIEW (time measurement of repaint methods)
#  -DDEBUG_UPDATE_ (output of update_(...) callers) 
-include Makefile_debug
