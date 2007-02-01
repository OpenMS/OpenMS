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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODELFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODELFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/KERNEL/DFeature.h>

namespace OpenMS
{

  /** @brief Abstract base class for ModelFitter-Module of FeatureFinder.
   
  		Every derived class has to implement the static functions
      "T* create()" and "const String getProductName()" (see FactoryProduct for details)
      
      @ingroup FeatureFinder
    
  */
  class BaseModelFitter 
    : public FeaFiModule
  {
  	public:
	    /** @brief Inner Classes for Exception handling
	     		
	     		UnableToFit-Excpetion if ModelFitter can not fit a model
					i.e. data set with standard deviation of zero 
			*/
	   	class UnableToFit
	     : public Exception::Base
			{
				public:
				 UnableToFit(const char* file, int line, const char* function, const std::string& name , const std::string& message) throw();
				 virtual ~UnableToFit() throw();
			};
	
	    /// Default constructor
	    BaseModelFitter();
	
	    /// copy constructor
	    BaseModelFitter(const BaseModelFitter& source);
	
	    /// destructor
	    virtual ~BaseModelFitter();
	
	    /// assignment operator
	    virtual BaseModelFitter& operator = (const BaseModelFitter& source);
	
	    /// register all derived classes here 
	    static void registerChildren();
	
	    /** @brief fit the range of indices onto a model
	    
	       @param extension range of peaks that ought to be fitted 
	       @return feature 
	    */
	    virtual DFeature<2> fit(const IndexSet& extension) throw(UnableToFit) =0;

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODELFITTER_H
