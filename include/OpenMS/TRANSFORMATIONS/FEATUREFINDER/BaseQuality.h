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
// $Id: BaseQuality.h,v 1.7 2006/04/25 12:41:44 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEQUALITY_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEQUALITY_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	class IndexSet;

  /** @brief Abstract base class to calculate quality for Modelfitting
    
  	 Abstract base class to calculate quality for Modelfitting
      every derived class has to implement the static functions
      "T* create()" and "const String getName()" (see FactoryProduct for details)
      
      @ingroup FeatureFinder
      
   */
  class BaseQuality: public FeaFiModule
  {
  public:
    /// standard constructor
    BaseQuality();

    /// copy constructor
    BaseQuality(const BaseQuality& source);

    /// destructor 
    virtual ~BaseQuality();

    /// assignment operator
    virtual BaseQuality& operator = (const BaseQuality& source);

    /// register all derived classes here 
    static void registerChildren();

    /** @brief calculate probability of peaks given by @p set given the model @p model
		
				This quality will be maximized by the ModelFitter 
		*/
    virtual double evaluate(const IndexSet& set, const BaseModel<2>& model)=0;
    
    /** @brief calculate probability of peaks given by @p set given the model @p model
        along dimension @p dim.
     */
    virtual double evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)=0;
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEQUALITY_H
