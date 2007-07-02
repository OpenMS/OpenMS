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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <queue>
#include <algorithm>
#include <cmath>

namespace OpenMS
{

	/**
	  @brief This extender does not extend the Seeding region.
		
		This module simply returns the seeding region without any further extension.
		It is used if your Seeding already returns an extended region, or for testing purposes.

		@ingroup FeatureFinder
	*/
  class DummyExtender 
    : public BaseExtender
  {

  public:
  	/// Default constructor
    DummyExtender();

    /// destructor
    virtual ~DummyExtender();

    /// Copy constructor
    DummyExtender(const DummyExtender& rhs);
    
    /// Assignment operator
    DummyExtender& operator= (const DummyExtender& rhs);

    /// return next seed
    const ChargedIndexSet& extend(const ChargedIndexSet& seed_region);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new DummyExtender();
    }

		/// returns the name of this module 
    static const String getProductName()
    {
      return "DummyExtender";
    }
                  
  protected:
		/// initialize member
  	virtual void updateMembers_();
  					
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYEXTENDER_H
