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
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include ../config_defs.mak
DIRECTORY=KERNEL

ADD_CXXFLAGS=-D_FILE_OFFSET_BITS=64

THISLIB=$(LIBNAME)

CPP_SOURCES= \
	KernelTraits.C\
	DSpectrum.C\
	DFeature.C\
	DimensionDescription.C\
	DPeak.C\
	DPickedPeak.C\
	DPeakArray.C\
	DRawDataPoint.C\
	StandardTypes.C\
	DFeatureMap.C\
	MSSpectrum.C\
	MSExperiment.C\
	MSExperimentExtern.C\
	DPeakConstReferenceArray.C\
	ConsensusFeature.C\
	ConsensusPeak.C\
	AreaIterator.C

include ../common.mak
-include .Dependencies
