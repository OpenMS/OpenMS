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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASESEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASESEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/DATASTRUCTURES/IndexSet.h>

namespace OpenMS
{
  /** @brief Abstract base class for Seeder-Module of FeatureFinder 
   
  	 	Abstract base class for Seeder-Module of FeatureFinder
      every derived class has to implement the static functions
      "T* create()" and "const String getName()" (see FactoryProduct for details)
      
      @ingroup FeatureFinder
      
   */
  class BaseSeeder: public FeaFiModule
  {

  public:
    /// standard constructor 
    BaseSeeder();

    /// copy constructor 
    BaseSeeder(const BaseSeeder& source);

    /// destructor 
    virtual ~BaseSeeder();

    /// assignment operator 
    virtual BaseSeeder& operator = (const BaseSeeder& source);

    /// register all derived classes here 
    static void registerChildren();

    /// return next seed 
    virtual IndexSet nextSeed() throw (NoSuccessor)=0;

    
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASESEEDER_H
