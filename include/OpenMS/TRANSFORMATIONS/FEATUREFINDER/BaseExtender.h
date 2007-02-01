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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/KERNEL/KernelTraits.h>

namespace OpenMS
{

  /**  
  	@brief Abstract base class for Extender-Module of FeatureFinder.
    
		Every derived class has to implement the static functions
    "T* create()" and "const String getProductName()" (see FactoryProduct for details).
    
    @ingroup FeatureFinder
  */
  class BaseExtender 
    : public FeaFiModule
  {

	  public:
	    /// Default constructor
	    BaseExtender();
	
	    /// copy constructor  
	    BaseExtender(const BaseExtender& source);
	
	    ///   destructor  
	    virtual ~BaseExtender();
	
	    /// assignment operator  
	    virtual BaseExtender& operator = (const BaseExtender& source);
	
	    /// register all derived classes here  
	    static void registerChildren();
	
	    /** 
	    	@brief extend given seed  
	    
				@param seed_region index of peak that serves as seed for feature
				@return IndexSet of peaks that could be part of a feature (only valid until next call to extend) 
	    */
	    virtual const IndexSet& extend(const IndexSet& seed_region)=0;
				
	  protected:
	    
	    /// Set of indizes representing the region that belongs to this feature.
	    IndexSet region_;
    
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEEXTENDER_H
