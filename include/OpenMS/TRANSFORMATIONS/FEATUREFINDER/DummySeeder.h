// -*- C++: make; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <vector>
#include <iostream>

namespace OpenMS
{

	/** @brief Dummy seeder class which always returns the same seed (1)
	
		This class can be used for Feature Finder instances that don't need a seeding module.
		One example is the class SweepExtender which performs Seeding and Extension
		at the same time.
	 	
		@ingroup FeatureFinder
		
	*/ 
  class DummySeeder 
    : public BaseSeeder
  {

  public:
	  	
    /// standard constructor
    DummySeeder();

    /// destructor 
    virtual ~DummySeeder();

    /// return next seed (always 1)
    Index nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new DummySeeder();
    }

    static const String getName()
    {
      return "DummySeeder";
    }

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H
