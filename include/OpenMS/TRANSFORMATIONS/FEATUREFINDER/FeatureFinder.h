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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
	/**
		@brief The base feature finding class
		
		This class calls a FeatureFinderAlgorithm
		
		@ingroup FeatureFinder
	*/
	class FeatureFinder
		: public ProgressLogger,
			public DefaultParamHandler,
			public FeatureFinderDefs
	{
		public:
			/// Default constructor.
	    FeatureFinder();
	    
	    /// destructor
	    virtual ~FeatureFinder();
			
			/// Execures the FeatureFinder using the given algorithm
			template<class PeakType, class FeatureType>
			void run(MSExperiment<PeakType> map, FeatureMap<FeatureType> features);
						
			/// Returns a non-mutable reference to a peak flag
	    inline const Flag& getPeakFlag(const IDX& index) const
	    {
	    	return flags_[index.first][index.second];
	    }
			
			/// Returns mutable reference to a peak flag
	    inline Flag& getPeakFlag(const IDX& index) 
	    { 
	    	return flags_[index.first][index.second];
	    }

		protected:
			/// Flags map that corresponds to input data
	    std::vector< std::vector<Flag> > flags_;

	}; //class
} //namespace

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
namespace OpenMS
{
	/// Execures the FeatureFinder using the given algorithm
	template<class PeakType, class FeatureType>
	void FeatureFinder::run(MSExperiment<PeakType> map, FeatureMap<FeatureType> features)
	{
		//@todo check for empty map, empty scans, MS>1 scans, if ranges are updated (Marc, Clemens, Marcel)
		
		//resize peak flag vector
		flags_.resize(map.size());
		for (UInt i=0; i<map.size(); ++i)
		{
  		flags_[i].assign(map[i].size(), UNUSED);
		}
		
		String algorithm_name = param_.getValue("algorithm");
		FeatureFinderAlgorithm<PeakType, FeatureType>* algorithm = Factory<FeatureFinderAlgorithm<PeakType, FeatureType> >::create(algorithm_name);
		
		algorithm->setData(map,features,*this);
		algorithm->run();

		delete(algorithm);
	}	
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
