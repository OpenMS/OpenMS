// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/CONCEPT/Factory.h>

// all from BaseModelFitter derived classes
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakFitter.h>

namespace OpenMS
{
	void BaseModelFitter::registerChildren()
	{
		Factory<BaseModelFitter>::registerProduct(SimpleModelFitter::getProductName(), &SimpleModelFitter::create);
		Factory<BaseModelFitter>::registerProduct(ExtendedModelFitter::getProductName(), &ExtendedModelFitter::create);
		Factory<BaseModelFitter>::registerProduct(PeakFitter::getProductName(), &PeakFitter::create);
	}	

	BaseModelFitter::BaseModelFitter() 
		: FeaFiModule()
	{	
	}

	BaseModelFitter::BaseModelFitter(const BaseModelFitter& source)
		: FeaFiModule(source) 
	{
	}

	BaseModelFitter::~BaseModelFitter()
	{
	}

  
	BaseModelFitter& BaseModelFitter::operator = (const BaseModelFitter& source)
	{
		if (&source == this) return *this;
		
		FeaFiModule::operator = (source);
		
		return *this;
	}

	BaseModelFitter::UnableToFit::UnableToFit(const char* file, int line, const char* function, const std::string& name , const std::string& message) throw()
		:	Base(file, line, function, name, message)
	{
	}

	BaseModelFitter::UnableToFit::~UnableToFit() throw()
	{	
	}

}

