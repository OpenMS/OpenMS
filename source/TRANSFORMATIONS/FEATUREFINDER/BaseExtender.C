// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/CONCEPT/Factory.h>

// all from BaseExtender derived classes
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SweepExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PPExtender.h>
#ifdef GSL_DEF
	#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/WaveletExtender.h>
#endif

namespace OpenMS
{
	void BaseExtender::registerChildren()
	{
		Factory<BaseExtender>::registerProduct(SimpleExtender::getName(), &SimpleExtender::create);
		Factory<BaseExtender>::registerProduct(SweepExtender::getName(), &SweepExtender::create);
		Factory<BaseExtender>::registerProduct(PPExtender::getName(), &PPExtender::create);
#ifdef GSL_DEF
		Factory<BaseExtender>::registerProduct(WaveletExtender::getName(), &WaveletExtender::create);
#endif
	}	

	BaseExtender::BaseExtender() : FeaFiModule() 
	{
	}


	BaseExtender::BaseExtender(const BaseExtender& source)
		: FeaFiModule(source) {}

	BaseExtender::~BaseExtender(){}
	
	BaseExtender& BaseExtender::operator = (const BaseExtender& source)
	{
		FeaFiModule::operator = (source);
		return *this;
	}
}


