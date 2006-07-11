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
# $Maintainer: Ole Schulz-Trieglaff $
# --------------------------------------------------------------------------


include ../../config_defs.mak
DIRECTORY=TRANSFORMATIONS/FEATUREFINDER

THISLIB=$(LIBNAME)


CPP_SOURCES=\
	FeaFiModule.C\
	BaseExtender.C\
	BaseSeeder.C\
	BaseModelFitter.C\
	BaseModel.C\
	ModelDescription.C\
	SimpleExtender.C\
	SimpleSeeder.C\
	LessSimplerSeeder.C\
	SimpleModelFitter.C\
	GaussModel.C\
	ProductModel.C\
	FeatureFinder.C\
	FastExtender.C\
	RobustSeeder.C\
	BaseQuality.C\
	Correlation.C\
	IsotopeModel.C\
	BiGaussModel.C\
	FeaFiTraits.C
	
include ../../common.mak
-include .Dependencies

# Uncomment this line to get the FeaFi debug output
ADD_CXXFLAGS=-DDEBUG_FEATUREFINDER
 
