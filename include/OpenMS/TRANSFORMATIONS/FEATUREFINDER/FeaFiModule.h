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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/CONCEPT/Types.h>
#include <set>

namespace OpenMS
{
	class FeaFiTraits;

  /** 
  	@brief Class to hold a module of the FeatureFinder algorithm 
		module accesses datastructures using BaseFeaFiTraits.
      
		@ingroup FeatureFinder
  */
	class FeaFiModule 
    : public FactoryProduct
  {	
  	public:
			/// Index in a MSExperiment
			typedef std::pair<UnsignedInt,UnsignedInt> IDX;
			/// Index set
			typedef std::set<IDX> IndexSet;
		
			/** 
				@brief Inner Classes for Exception handling
			  
				NoSuccessor-Excpetion if getNext*** or getPrev***-Methods are called on an index 
				that has no successor or predecessor 
			*/
			class NoSuccessor
			 : public Exception::Base
			{
			public:
			 NoSuccessor(const char* file, int line, const char* function, const IDX& index) throw();
			 
			 virtual ~NoSuccessor() throw();
			 
			protected:
			 IDX index_;  // index without successor/predecessor
			};

			/// Default constructor 
			FeaFiModule();
			
			/// copy constructor 
			FeaFiModule(const FeaFiModule& source);
			
			/// destructor 
			virtual ~FeaFiModule();
			
			/// assignment operator 
			virtual FeaFiModule& operator = (const FeaFiModule& source);
			
			/// set FeatureFinder traits 
			void setTraits(FeaFiTraits* traits);

	  protected:
	  	/// Pointer to the tratis class
	   	FeaFiTraits* traits_;
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
