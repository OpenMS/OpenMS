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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>

namespace OpenMS
{
	template <typename,typename> class FeatureFinderAlgorithm;
	/**
		@brief The base class of the Feature Finding algorithm.
		
		@ingroup FeatureFinder
	*/
	template <class PeakType, class FeatureType>
	class FeatureFinder
		: public ProgressLogger,
			public DefaultParamHandler,
			public FeatureFinderDefs
	{
		friend class FeatureFinderAlgorithm<PeakType, FeatureType>;
		
		public:
			/// Input map type
			typedef MSExperiment<PeakType> MapType;
			/// Coordinate/Position type of peaks
			typedef typename PeakType::CoordinateType CoordinateType;
			/// Intensity type of peaks
			typedef typename PeakType::IntensityType IntensityType;
			/// Output feature type
			typedef FeatureMap<FeatureType> FeatureMapType;


			/// Default constructor.
	    FeatureFinder()
				: DefaultParamHandler("FeatureFinder"),
					map_(0),
					features_(0),
					flags_()
			{
			}
	    /// destructor
	    virtual ~FeatureFinder()
			{
			}
			/// Sets the input peak map
	    void setInput(const MapType& map)
		  {
				//TODO check for empty scans and MS!=1 scans	
	    	map_ = &map;
	    }
			/// Sets the output feature map
			void setOutput(FeatureMapType& features)
			{
				features_ = &features;
			}

			/// Execures the FeatureFinder using the given algorithm
	    void run();

		protected:
			/// Returns a non-mutable reference to the output feature map
			inline const FeatureMapType& getFeatureMap_() const 
			{ 
				return *features_; 
			}
			/// Returns a mutable reference to the output feature map
			inline FeatureMapType& getFeatureMap_() 
			{ 
				return *features_; 
			}
			/// Returns a non-mutable refernece to the input data
			inline const MapType& getData_() const 
			{ 
				return *map_; 
			}
			/// Returns a non-mutable reference to a peak flag
	    inline const Flag& getPeakFlag_(const IDX& index) const
	    {
	    	return flags_[index.first][index.second];
	    }
			/// Returns mutable reference to a peak flag
	    inline Flag& getPeakFlag_(const IDX& index) 
	    { 
	    	return flags_[index.first][index.second];
	    }
			/// Returns the intensity of a peak
	    inline IntensityType getPeakIntensity_(const IDX& index) const
	    { 
				//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_->size(), "Scan index outside of map->");
		    OPENMS_PRECONDITION(index.second<(*map_)[index.first].size(), "Peak index outside of scan!");
			
	    	return (*map_)[index.first][index.second].getIntensity(); 
	    }
	    /// Returns the  m/z of a peak
	    inline CoordinateType getPeakMz_(const IDX& index) const
	    { 
				//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_->size(), "Scan index outside of map->");
		    OPENMS_PRECONDITION(index.second<(*map_)[index.first].size(), "Peak index outside of scan!");
			
	    	return (*map_)[index.first][index.second].getMZ(); 
	    }
	    /// Returns the retention time of a peak
	    inline CoordinateType getPeakRt_(const IDX& index) const
	    { 
				//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_->size(), "Scan index outside of map->");
		    OPENMS_PRECONDITION(index.second<(*map_)[index.first].size(), "Peak index outside of scan!");
			
	    	return (*map_)[index.first].getRT();
	    }
	
	    /// fills @p index with the index of next peak in m/z dimension
	    inline void getNextMz_(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
		  	//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_->size(), "Scan index outside of map->");
		    OPENMS_PRECONDITION(index.second<(*map_)[index.first].size(), "Peak index outside of scan!");
    
	    	//At the last peak of this spectrum
	      if (index.second==(*map_)[index.first].size()-1)
	      {
	      	throw NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getNextMz", index);
	      }
	
	      ++index.second;
	    }
	
	    /// fills @p index with the index of previous peak in m/z dimension
	    inline void getPrevMz_(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
		  	//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_->size(), "Scan index outside of map->");
		    OPENMS_PRECONDITION(index.second<(*map_)[index.first].size(), "Peak index outside of scan!");
    
	      //begin of scan
	      if (index.second==0)
	      {
	      	throw NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getPrevMz", index);
				}
				
	      --index.second;
	    }
	


	/// fills @p index with the index of the nearest peak in the next scan
  void getNextRt_(IDX& index) throw (NoSuccessor, Exception::Precondition)
  {
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<map_->size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<(*map_)[index.first].size(), "Peak index outside of scan!");
		
		//last scan
    if (index.first == map_.size()-1)
    {
    	throw NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getNextRt", index);
		}
		
		// perform binary search to find the neighbour in rt dimension
		CoordinateType mz_pos = (*map_)[index.first][index.second].getMZ();	// mz value we want to find
		++index.first;
		typename MapType::SpectrumType::ConstIterator it = lower_bound((*map_)[index.first].begin(), (*map_)[index.first].end(), (*map_)[index.first-1][index.second], MapType::SpectrumType::PeakType::PositionLess());	
		
		// if the found peak is at the end of the spectrum, there is not much we can do...
		if ( it == (*map_)[index.first].end() )
		{
			// check for empty scans
			if ( (*map_)[index.first].size() > 0 )
	 			index.second = (*map_)[index.first].size()-1;
			else
				index.second = 0;
		}
		// if the found peak is at the beginning of the spectrum, there is also not much we can do ! 
		else if ( it == (*map_)[index.first].begin() ) 
		{
			index.second = 0;
		}
		// see if the next smaller one fits better
		else 
		{	
			// peak to the right is closer (in m/z dimension)
			if (it->getMZ() - mz_pos < mz_pos - (it-1)->getMZ() )
			{				
				index.second = it - (*map_)[index.first].begin(); 
			}
			else	// left one is closer
			{
				index.second = --it - (*map_)[index.first].begin(); 
			}
		}
  }
  
  
	/// fills @p index with the index of the nearest peak in the previous scan
	void getPrevRt_(IDX& index) throw (NoSuccessor, Exception::Precondition)
  {
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<(*map_)[index.first].size(), "Peak index outside of scan!");
		
		if (index.first>=map_.size() )
		{
			std::cout << "Scan index outside of map!" << std::endl;
			std::cout << index.first << " " << index.second << std::endl;
			return;
		}
		if (index.second>=(*map_)[index.first].size())
		{
			std::cout << "Peak index outside of scan!" << std::endl;
			std::cout << index.first << " " << index.second << std::endl;
			return;		
		}	
		
		// first scan
		if (index.first == 0)
		{
			throw NoSuccessor(__FILE__, __LINE__, "FeatureFinder::getPrevRt", index);
		}
		
		// perform binary search to find the neighbour in rt dimension
		CoordinateType mz_pos = (*map_)[index.first][index.second].getMZ();
		--index.first;
		typename MapType::SpectrumType::ConstIterator it = lower_bound((*map_)[index.first].begin(), 
		                                                                                  (*map_)[index.first].end(), 
																																								      (*map_)[index.first+1][index.second], 
																																								      MapType::SpectrumType::PeakType::PositionLess());	
		
		// if the found peak is at the end of the spectrum, there is not much we can do.
		if ( it == (*map_)[index.first].end() )
		{
	 		// check for empty scans
			if ( (*map_)[index.first].size() > 0 )
	 			index.second = (*map_)[index.first].size()-1;
			else
				index.second = 0;
		}
		// if the found peak is at the beginning of the spectrum, there is not much we can do.
		else if ( it == (*map_)[index.first].begin() ) 
		{
			index.second = 0;
		}
		// see if the next smaller one fits better
		else 
		{	
			// peak to the right is closer (in m/z dimension)
			if (it->getMZ() - mz_pos < mz_pos - (it-1)->getMZ() )
			{
				index.second = it - (*map_)[index.first].begin(); 
			}
			else
			{
				index.second = --it - (*map_)[index.first].begin(); 
			}
		}
  }
	
	///Calculates the convex hull of a index set and adds it to the feature
	void addConvexHull_(const IndexSet& set, Feature& f) const
	{
		std::vector< DPosition<2> > points;
		points.reserve(set.size());
		DPosition<2> tmp;
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
    {
    	tmp[RawDataPoint2D::MZ] = (*map_)[it->first][it->second].getMZ();
    	tmp[RawDataPoint2D::RT] = (*map_)[it->first].getRT();
    	points.push_back(tmp);
    }
		f.getConvexHulls().resize(f.getConvexHulls().size()+1);
		// assignment operator computes convex hull
		f.getConvexHulls()[f.getConvexHulls().size()-1] = points;	
	}

		protected:
			const MSExperiment<PeakType>* map_;
			FeatureMap<FeatureType>* features_;
	    std::vector< std::vector<Flag> > flags_;
	}; //class
} //namespace


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

namespace OpenMS
{
	template<class PeakType, class FeatureType>
	void FeatureFinder<PeakType,FeatureType>::run()
	{
		String algorithm_name = param_.getValue("algorithm");
		FeatureFinderAlgorithm<PeakType, FeatureType>* algorithm = Factory<FeatureFinderAlgorithm<PeakType, FeatureType> >::create(algorithm_name);
				
		algorithm->setFeatureFinder(*this);			
		algorithm->run();

		delete(algorithm);
	}
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
