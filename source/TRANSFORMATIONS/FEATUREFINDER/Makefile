# -*- Mode: makefile; tab-width: 2; -*-
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
# $Maintainer: Ole Schulz-Trieglaff, Clemens Groepl $
# --------------------------------------------------------------------------


DIRECTORY=TRANSFORMATIONS/FEATUREFINDER

THISLIB=$(LIBNAME)

include ../../config_defs.mak

CPP_SOURCES=\
	BaseExtender.C\
	BaseModel.C\
	BaseModelFitter.C\
	BaseQuality.C\
	BaseSeeder.C\
	BiGaussModel.C\
	Correlation.C\
	DummySeeder.C\
	FeaFiModule.C\
	FeaFiTraits.C\
	FeatureFinder.C\
	GaussModel.C\
	IsotopeModel.C\
	ModelDescription.C\
	ProductModel.C\
	SimpleExtender.C\
	SimpleModelFitter.C\
	SimpleSeeder.C\
	SweepExtender.C\
	PPExtender.C\
	WaveletSeeder.C\
	WaveletExtender.C

ifeq ($(GSL_SUPPORT),true)
 GSL_SOURCES=\
	EmgModel.C\
	ExtendedModelFitter.C\
	LmaGaussModel.C\
	LogNormalModel.C
endif

ADD_CPP_SOURCES = ${GSL_SOURCES} 

include ../../common.mak
-include .Dependencies

# Uncomment this line to get the FeaFi debug output
#ADD_CXXFLAGS=-DDEBUG_FEATUREFINDER
