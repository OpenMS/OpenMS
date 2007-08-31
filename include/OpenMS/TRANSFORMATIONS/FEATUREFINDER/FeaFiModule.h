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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

namespace OpenMS
{
	class FeaFiTraits;

  /** 
  	@brief Implements a module of the FeatureFinder algorithm.
      
		@ingroup FeatureFinder
  */
	template<class PeakType, class FeatureType>
	class FeaFiModule 
    : public DefaultParamHandler
  {	
  	public:
			/// FeatureFinder algorithm type
			typedef FeatureFinder<PeakType, FeatureType> FeatureFinderType;
			/// Input map type
			typedef typename FeatureFinderType::MapType MapType;
			/// Output Feature type
			typedef typename FeatureFinderType::FeatureMapType FeatureMapType;
			/// Coordinate/Position type
			typedef typename FeatureFinderType::CoordinateType CoordinateType;
			/// Intensity type
			typedef typename FeatureFinderType::IntensityType IntensityType;

			/// Default constructor 
			FeaFiModule()
			: DefaultParamHandler("FeaFiModule"), 
				traits_(0)
			{
			}
			/// copy constructor 
			FeaFiModule(const FeaFiModule& source)
			: DefaultParamHandler(source),
				traits_(source.traits_)
			{
			}
			
			/// destructor 
			virtual ~FeaFiModule()
			{
			}
			
			/// assignment operator 
			virtual FeaFiModule& operator = (const FeaFiModule& source)
			{
				if (&source == this) return *this;
		
				DefaultParamHandler::operator = (source);
				traits_ = source.traits_;
		
				return *this;
			}
			
			/// set FeatureFinder algorithm
			void setAlgorithm(FeatureFinderAlgorithm<PeakType,FeatureType>* traits)
			{
				traits_ = traits;
			}

	  protected:
	  	/// Pointer to the tratis class
	   	FeatureFinderAlgorithm<PeakType,FeatureType>* traits_;
	};
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
