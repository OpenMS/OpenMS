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
# $Maintainer:  $
# --------------------------------------------------------------------------


EXECUTABLE=SpecView_SpecAnnotate
DIRECTORY=APPLICATIONS/SPECVIEW

#default: all

all: $(EXECUTABLE)

#replace library-specific flags by application-
#specific compiler flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS) 


include ../../common.mak
-include .Dependencies

clean: default_clean
	@-$(RM) -rf .Dependencies
	@-$(RM) -rf ./SPECANNOTATE/FUNCTION/*~
	@-$(RM) -rf ./SPECANNOTATE/VISUAL/*~
	@-$(RM) -rf ./SPECANNOTATE/VISUAL/XML/*~
	@-$(RM) -rf ./SPECANNOTATE/VISUAL/DATABASE/*~
	@-$(RM) -rf ./SPECANNOTATE/*~

ADD_CXXFLAGS= -DQT_THREAD_SUPPORT -DQT_SHARED -D_REENTRANT -pipe -Wall -W -O2 -Wno-deprecated -DQT_NO_DEBUG 

OPENMS_INCLUDES=$(GUI_INCLUDES) -ISPECANNOTATE/FUNCTION -ISPECANNOTATE/VISUAL -ISPECANNOTATE/VISUAL/DATABASE 
LIBS=$(GUI_LIBS)

CPP_SOURCES=SPECANNOTATE/VISUAL/SpecAnnotate.C \
	SPECANNOTATE/VISUAL/AnnotateThread.C \
	SPECANNOTATE/FUNCTION/MySQLAdapter.C \
	SPECANNOTATE/FUNCTION/AminoAcid.C \
	SPECANNOTATE/FUNCTION/Enzyme.C \
	SPECANNOTATE/FUNCTION/Modification.C \
	SPECANNOTATE/FUNCTION/ProteinDigest.C \
	SPECANNOTATE/FUNCTION/Sample.C \
	SPECANNOTATE/FUNCTION/ModificationStringParser.C \
	SPECANNOTATE/FUNCTION/Annotation.C \
	SPECANNOTATE/VISUAL/SpectrumMDIWindowEnhanced.C \
	SPECANNOTATE/VISUAL/SampleDialog.C \
	SPECANNOTATE/VISUAL/SettingsDialog.C \
	SPECANNOTATE/VISUAL/Annotate.C \
	SPECANNOTATE/VISUAL/InputModifications.C \
	SPECANNOTATE/VISUAL/DATABASE/EditModification.C \
	SPECANNOTATE/VISUAL/DATABASE/EditEnzyme.C \
	SPECANNOTATE/VISUAL/DATABASE/EditAminoacid.C \
	SPECANNOTATE/VISUAL/DATABASE/EditProtein.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewRealizedModificationPositionless.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewSequence.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewDigestFragment.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewModificationCombination.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewModificationCombinationPositionless.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewProteinModificationScenario.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewRealizedModification.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewSample.C \
	SPECANNOTATE/VISUAL/DATABASE/ViewAnnotation.C \
	SPECANNOTATE/VISUAL/XML/EditEnzymeXML.C \
	SPECANNOTATE/VISUAL/XML/EditModificationXML.C \
	SpecView_SpecAnnotate.C

MOC_SOURCES= SPECANNOTATE/VISUAL/moc_SpecAnnotate.C \
	SPECANNOTATE/VISUAL/moc_SpectrumMDIWindowEnhanced.C \
	SPECANNOTATE/VISUAL/moc_SampleDialog.C \
	SPECANNOTATE/VISUAL/moc_SettingsDialog.C \
	SPECANNOTATE/VISUAL/moc_Annotate.C \
	SPECANNOTATE/VISUAL/moc_InputModifications.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_EditModification.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_EditEnzyme.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_EditAminoacid.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_EditProtein.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewRealizedModificationPositionless.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewSequence.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewDigestFragment.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewModificationCombination.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewModificationCombinationPositionless.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewProteinModificationScenario.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewRealizedModification.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewSample.C \
	SPECANNOTATE/VISUAL/DATABASE/moc_ViewAnnotation.C \
	SPECANNOTATE/VISUAL/XML/moc_EditEnzymeXML.C \
	SPECANNOTATE/VISUAL/XML/moc_EditModificationXML.C 


SpecView_SpecAnnotate: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LIBS)


SPECANNOTATE/VISUAL/moc_SpecAnnotate.C: $(MOC) SPECANNOTATE/VISUAL/SpecAnnotate.h
	$(MOC) SPECANNOTATE/VISUAL/SpecAnnotate.h -o SPECANNOTATE/VISUAL/moc_SpecAnnotate.C

SPECANNOTATE/VISUAL/moc_SpectrumMDIWindowEnhanced.C: $(MOC) SPECANNOTATE/VISUAL/SpectrumMDIWindowEnhanced.h
	$(MOC) SPECANNOTATE/VISUAL/SpectrumMDIWindowEnhanced.h -o SPECANNOTATE/VISUAL/moc_SpectrumMDIWindowEnhanced.C

SPECANNOTATE/VISUAL/moc_SampleDialog.C: $(MOC) SPECANNOTATE/VISUAL/SampleDialog.h
	$(MOC) SPECANNOTATE/VISUAL/SampleDialog.h -o SPECANNOTATE/VISUAL/moc_SampleDialog.C

SPECANNOTATE/VISUAL/moc_SettingsDialog.C: $(MOC) SPECANNOTATE/VISUAL/SettingsDialog.h
	$(MOC) SPECANNOTATE/VISUAL/SettingsDialog.h -o SPECANNOTATE/VISUAL/moc_SettingsDialog.C

SPECANNOTATE/VISUAL/moc_Annotate.C: $(MOC) SPECANNOTATE/VISUAL/Annotate.h
	$(MOC) SPECANNOTATE/VISUAL/Annotate.h -o SPECANNOTATE/VISUAL/moc_Annotate.C

SPECANNOTATE/VISUAL/moc_InputModifications.C: $(MOC) SPECANNOTATE/VISUAL/InputModifications.h
	$(MOC) SPECANNOTATE/VISUAL/InputModifications.h -o SPECANNOTATE/VISUAL/moc_InputModifications.C

SPECANNOTATE/VISUAL/DATABASE/moc_EditModification.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/EditModification.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/EditModification.h -o SPECANNOTATE/VISUAL/DATABASE/moc_EditModification.C

SPECANNOTATE/VISUAL/DATABASE/moc_EditEnzyme.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/EditEnzyme.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/EditEnzyme.h -o SPECANNOTATE/VISUAL/DATABASE/moc_EditEnzyme.C

SPECANNOTATE/VISUAL/DATABASE/moc_EditAminoacid.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/EditAminoacid.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/EditAminoacid.h -o SPECANNOTATE/VISUAL/DATABASE/moc_EditAminoacid.C

SPECANNOTATE/VISUAL/DATABASE/moc_EditProtein.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/EditProtein.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/EditProtein.h -o SPECANNOTATE/VISUAL/DATABASE/moc_EditProtein.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewRealizedModificationPositionless.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewRealizedModificationPositionless.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewRealizedModificationPositionless.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewRealizedModificationPositionless.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewSequence.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewSequence.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewSequence.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewSequence.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewDigestFragment.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewDigestFragment.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewDigestFragment.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewDigestFragment.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewModificationCombination.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewModificationCombination.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewModificationCombination.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewModificationCombination.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewModificationCombinationPositionless.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewModificationCombinationPositionless.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewModificationCombinationPositionless.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewModificationCombinationPositionless.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewProteinModificationScenario.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewProteinModificationScenario.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewProteinModificationScenario.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewProteinModificationScenario.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewRealizedModification.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewRealizedModification.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewRealizedModification.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewRealizedModification.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewSample.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewSample.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewSample.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewSample.C

SPECANNOTATE/VISUAL/DATABASE/moc_ViewAnnotation.C: $(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewAnnotation.h
	$(MOC) SPECANNOTATE/VISUAL/DATABASE/ViewAnnotation.h -o SPECANNOTATE/VISUAL/DATABASE/moc_ViewAnnotation.C

SPECANNOTATE/VISUAL/XML/moc_EditEnzymeXML.C: $(MOC) SPECANNOTATE/VISUAL/XML/EditEnzymeXML.h
	$(MOC) SPECANNOTATE/VISUAL/XML/EditEnzymeXML.h -o SPECANNOTATE/VISUAL/XML/moc_EditEnzymeXML.C

SPECANNOTATE/VISUAL/XML/moc_EditModificationXML.C: $(MOC) SPECANNOTATE/VISUAL/XML/EditModificationXML.h
	$(MOC) SPECANNOTATE/VISUAL/XML/EditModificationXML.h -o SPECANNOTATE/VISUAL/XML/moc_EditModificationXML.C
