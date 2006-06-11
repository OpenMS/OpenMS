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
// $Id: BaseQuality.C,v 1.4 2006/03/28 10:07:05 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff  $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseQuality.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/DATASTRUCTURES/IndexSet.h>

// all from BaseQuality derived classes
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Correlation.h>

namespace OpenMS
{
	void BaseQuality::registerChildren()
	{
		Factory<BaseQuality>::registerProduct("Correlation", &Correlation::create);
	}

	BaseQuality::BaseQuality(): FeaFiModule(){}

	BaseQuality::BaseQuality(const BaseQuality& source)
		: FeaFiModule(source)
	{}

	BaseQuality::~BaseQuality() {}

  BaseQuality& BaseQuality::operator = (const BaseQuality& source)
	{
		FeaFiModule::operator = (source);
		return *this;
	}
}
