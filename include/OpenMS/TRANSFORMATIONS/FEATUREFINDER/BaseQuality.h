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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEQUALITY_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEQUALITY_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  /** @brief Abstract base class to calculate quality for Modelfitting
    
  	 Abstract base class to calculate quality for Modelfitting
      every derived class has to implement the static functions
      "T* create()" and "const String getProductName()" (see FactoryProduct for details)
      
      @ingroup FeatureFinder
      
   */
  class BaseQuality: public FeaFiModule
  {
	  public:
	    /// Default constructor
	    BaseQuality();
	
	    /// copy constructor
	    BaseQuality(const BaseQuality& source);
	
	    /// destructor 
	    virtual ~BaseQuality();
	
	    /// assignment operator
	    virtual BaseQuality& operator = (const BaseQuality& source);
	
	    /// register all derived classes here 
	    static void registerChildren();
	
	    /// Computes the goodness of fit between model @p model and peaks denoted by @p set
	    virtual double evaluate(const IndexSet& set, const BaseModel<2>& model)=0;
	    
	    /// Computes the goodness of fit between model @p model and peaks @p set along dimension @p dim.
	    virtual double evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)=0;
			
			/// Returns the significance of the last fit (if applicable, otherwise -1
			double getPvalue()
			{ 
				return pval_;
			}
			
		protected:
			/// Significance of fit
			double pval_;
			
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEQUALITY_H
