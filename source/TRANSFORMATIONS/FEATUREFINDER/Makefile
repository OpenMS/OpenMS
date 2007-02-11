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
# $Maintainer: Ole Schulz-Trieglaff, Clemens Groepl $
# --------------------------------------------------------------------------

DIRECTORY=TRANSFORMATIONS/FEATUREFINDER

THISLIB=$(LIBNAME)

include ../../config_defs.mak

CPP_SOURCES=\
	FeaFiModule.C\
	ModelDescription.C\
	BaseModel.C\
	ProductModel.C\
	GaussModel.C\
	BiGaussModel.C\
	IsotopeModel.C\
	EmgModel.C\
	LmaGaussModel.C\
	LogNormalModel.C\
	BaseExtender.C\
	SimpleExtender.C\
	BaseModelFitter.C\
	SimpleModelFitter.C\
	ExtendedModelFitter.C\
	BaseSeeder.C\
	SimpleSeeder.C\
	PickedPeakSeeder.C\
	MarrWaveletSeeder.C\
	IsotopeWaveletSeeder.C\
	BaseQuality.C\
	Correlation.C\
	RankCorrelation.C\
	EuclideanDistance.C\
	FeaFiTraits.C\
	FeatureFinder.C\
	PeakSeeder.C\
	PeakExtender.C\
	PeakFitter.C

ADD_CPP_SOURCES = $(CGAL_SOURCES)

include ../../common.mak
-include .Dependencies

# For debugging the FeatureFinder, please create the file 'Makefile_debug'
# and add the line 'ADD_CXXFLAGS=-DDEBUG_FEATUREFINDER' to it
-include Makefile_debug
