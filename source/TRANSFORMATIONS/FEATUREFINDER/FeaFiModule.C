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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>

namespace OpenMS
{
	FeaFiModule::FeaFiModule()
		: FactoryProduct("FeaFiModule"), 
			traits_(0)
	{
	}

	FeaFiModule::FeaFiModule(const FeaFiModule& source)
		:FactoryProduct(source),
		traits_(source.traits_)
	{
	}

	FeaFiModule::~FeaFiModule()
	{
	}

	FeaFiModule& FeaFiModule::operator = (const FeaFiModule& source)
	{
		if (&source == this) return *this;
		
		FactoryProduct::operator = (source);
		traits_ = source.traits_;
		
		return *this;
	}

	void FeaFiModule::setTraits(FeaFiTraits* traits)
	{
		traits_ = traits;
	}

	FeaFiModule::NoSuccessor::NoSuccessor(const char* file, int line, const char* function, const IDX& index) throw()
		:	Base(file, line, function, "NoSuccessor", "no successor/predecessor"), 
			index_(index)
	{
		what_ = String("there is no successor/predecessor for the given Index: ") + index_.first + "/" + index_.second;
		Exception::globalHandler.setMessage(what_);
	}

	FeaFiModule::NoSuccessor::~NoSuccessor() throw()
	{
	}
}
