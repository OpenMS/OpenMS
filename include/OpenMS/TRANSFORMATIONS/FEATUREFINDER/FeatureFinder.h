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
			public FeatureFinderDefs
	{
		public:
			/// Default constructor.
	    FeatureFinder();
	    
	    /// destructor
	    virtual ~FeatureFinder();
			
			/**
				@brief Executes the FeatureFinder using the given algorithm
				
				Note that there are several constraints for the @p map. It must only contain MS1 scans and
				you have to call updateRanges() before passing it to this method.
				
				@param algorithm_name Name of the feature finding algorithm to use
				@param map Input peak map
				@param features Output feature map
				@param param Algorithm parameters
			*/
			template<class PeakType, class FeatureType>
			void run(const String& algorithm_name, MSExperiment<PeakType> const & map, FeatureMap<FeatureType> & features, const Param& param);
						
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
			
			/// Returns the default parameters for a the algorithm with name @p algorithm_name
			Param getParameters(const String& algorithm_name) const;
			
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
	void FeatureFinder::run(const String& algorithm_name, MSExperiment<PeakType> const & map, FeatureMap<FeatureType> & features, const Param& param)
	{
		//Nothing to do if there is no data
		if (map.size()==0)
		{
			return;
		}
		//We need updated ranges => check number of peaks
		else if (map.getSize()==0)
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"FeatureFinder needs updated ranges on input map. Aborting!");
		}
		//We need MS1 data only => check levels
		else if (map.getMSLevels().size()!=1 || map.getMSLevels()[0]!=1)
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"FeatureFinder can only operate on MS level 1 data. Please do not use MS/MS data. Aborting!");
		}
		//@todo Should we check for empty scans? I think the FF should be able to deal with them! (Marc, Clemens, Marcel)
		
		//resize peak flag vector
		flags_.resize(map.size());
		for (UInt i=0; i<map.size(); ++i)
		{
  		flags_[i].assign(map[i].size(), UNUSED);
		}
		
		if (algorithm_name!="none")
		{
			FeatureFinderAlgorithm<PeakType, FeatureType>* algorithm =
				Factory<FeatureFinderAlgorithm<PeakType, FeatureType> >::create(algorithm_name);
			algorithm->setParameters(param);
			algorithm->setData(map,features,*this);
			algorithm->run();
			delete(algorithm);
		}	
	}
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
