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

DIRECTORY=EXAMPLES

CPP_SOURCES=\
  Tutorial_DPosition.C\
  Tutorial_DRange.C\
	Tutorial_MSSpectrum.C\
	Tutorial_MSExperiment.C\
	Tutorial_FeatureMap.C\
	Tutorial_RangeManager.C\
	Tutorial_MetaInfo.C\
	Tutorial_FileIO.C\
	Tutorial_Param.C\
	Tutorial_PeakPickerCWT.C\
	Tutorial_TopHatFilter.C\
	Tutorial_GaussFilter.C\
	Tutorial_SavitzkyGolayFilter.C\
	Tutorial_FeatureFinder.C\
	Tutorial_MapAlignment.C\
	Tutorial_InternalCalibration.C\
	Tutorial_TOFCalibration.C\
	Tutorial_Residue.C\
	Tutorial_Element.C\
	Tutorial_EmpiricalFormula.C\
	Tutorial_AASequence.C\
	Tutorial_TheoreticalSpectrumGenerator.C\
	Tutorial_Spectrum1D.C\
	Tutorial_ParamEditor.C\
	Tutorial_Unlabeled.C\
	Tutorial_Labeled.C\
	Tutorial_Clustering.C\
	Tutorial_typeAsString.C

OBJECTS=$(SOURCES:.C=.o)

help:
	@echo ""
	@echo "Make targets are:"
	@echo "all -- builds all executables"
	@echo "depend -- updates dependencies"
	@echo "clean -- cleans up the directory"
	@echo "execute -- executes all examples"
	@echo ""
	@echo "Internal targets are:"
	@echo "ex -- Builds the executable 'example' from 'example.C'"
	@echo ""

include ../common.mak
-include .Dependencies

.o:
	$(CXX) $(CXXFLAGS) $*.o $(LDFLAGS) -o $@ $(LIBS)

EXECUTABLES=$(SOURCES:.C=)

all:	$(EXECUTABLES)

clean:	default_clean
	@rm -rf output/Tutorial_*

execute: $(EXECUTABLES)
	@for i in $(EXECUTABLES); do\
		echo ;\
		echo "===================================================================================";\
		echo "./$$i";\
		./$$i;\
	done

#dummy target for little test programs
example: example.o

#Main example target
#as example is not in the 'depend' and 'clean' target, we need to clean up before building it
ex:
	@rm -rf example.o example
	@make example

dummy:
