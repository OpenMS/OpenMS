// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FITTER1D_H

#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>

namespace OpenMS
{
	class InterpolationModel;

	/**
	@brief Abstract base class for all 1D-dimensional model fitter.

	Every derived class has to implement the static functions
	"T* create()" and "const String getProductName()" (see DefaultParamHandler for details)
		
	@htmlinclude OpenMS_Fitter1D.parameters
  
	@ingroup FeatureFinder
	*/
	class OPENMS_DLLAPI Fitter1D
    : public DefaultParamHandler,
      public FeatureFinderDefs
	{
	 public:
            
		/// IndexSet
		typedef IsotopeCluster::IndexSet IndexSet;
		/// IndexSet with charge information
		typedef IsotopeCluster::ChargedIndexSet ChargedIndexSet;
		/// Single coordinate
		typedef Feature::CoordinateType CoordinateType;
		/// Quality of a feature
		typedef Feature::QualityType QualityType;
		/// Raw data point type
		typedef Peak1D PeakType;
		/// Raw data container type using for the temporary storage of the input data
		typedef std::vector<PeakType > RawDataArrayType;
		/// Raw data iterator
		typedef RawDataArrayType::iterator PeakIterator;
          
		/// Default constructor.
		Fitter1D();
            
		/// copy constructor
		Fitter1D(const Fitter1D& source);
    
		/// destructor
		virtual ~Fitter1D()
		{
		}
    
		/// assignment operator
		virtual Fitter1D& operator = (const Fitter1D& source);
                            
		/// return interpolation model
		virtual QualityType fit1d(const RawDataArrayType& /* range */, InterpolationModel*& /* model */)
		{
			throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
    
		/// register all derived classes here
		static void registerChildren();
		
	 protected:
		
		/// standard derivation in bounding box		
		CoordinateType tolerance_stdev_box_;
		/// minimum of the bounding box
		CoordinateType min_;
		/// maximum of the bounding box
		CoordinateType max_;
		/// standard derivation
		CoordinateType stdev1_; 
		/// standard derivation            
		CoordinateType stdev2_;
		/// basic statistics           
		Math::BasicStatistics<> statistics_;
		/// interpolation step size
		CoordinateType interpolation_step_;
            
		virtual void updateMembers_();
				
	};

}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FITTER1D_H
