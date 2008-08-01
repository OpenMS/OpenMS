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
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include ../config_defs.mak
DIRECTORY=KERNEL

ADD_CXXFLAGS=-D_FILE_OFFSET_BITS=64

CPP_SOURCES= \
	PeakIndex.C\
  DPeak.C\
  Peak1D.C\
  Peak2D.C\
	DRichPeak.C\
  RichPeak1D.C\
  RichPeak2D.C\
	Feature.C\
  DPeakArray.C\
  StandardTypes.C\
	RangeManager.C\
	FeatureMap.C\
	DSpectrum.C\
  MSSpectrum.C\
  AreaIterator.C\
	MSExperiment.C\
	MSExperimentExtern.C\
	FeatureHandle.C\
	ConsensusFeature.C\
	ConsensusMap.C

include ../common.mak
-include .Dependencies

