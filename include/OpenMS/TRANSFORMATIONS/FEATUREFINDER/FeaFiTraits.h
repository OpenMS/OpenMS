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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>

#include <set>

namespace OpenMS
{
	class BaseSeeder;
	class BaseExtender;
	class BaseModelFitter;
	/**
		@brief Traits class for the feature finding algorithm.
		
		This class is rather an "umbrella" for the different modules / steps of the algorithm
		than a traits class in the traditional sense.
		
		@todo Allow choice between different map types (MSExperiment or MSExperimentExtern) (Ole)
		
		@ingroup FeatureFinder 	
	*/
	class FeaFiTraits
	{
		public:
			/// Index in a MSExperiment ( first index denotes rt, second m/z )
			typedef std::pair<UnsignedInt,UnsignedInt> IDX;
			/// Index set
			typedef std::set<IDX> IndexSet;
	
	    /// Defines the coordinates of peaks / features.
	    enum DimensionId
	    {
        RT = DimensionDescription < LCMS_Tag >::RT,
        MZ = DimensionDescription < LCMS_Tag >::MZ
	    };
	
	    /// Flag for each data point
	    enum Flag { UNUSED, SEED, INSIDE_FEATURE };
			
			/// Internal map type
	    typedef MSExperimentExtern< DPeak<1> > MapType;
			/// Intensity type of the map
	    typedef MapType::IntensityType IntensityType;
	    /// Coordinate type of the map
	    typedef MapType::CoordinateType CoordinateType;
				
	    /// 2D position type (needed for models)
	    typedef DPosition<2> PositionType2D;
	
			/// No successor exception
	    typedef FeaFiModule::NoSuccessor NoSuccessor;
	
	    /// Default constructor
	    FeaFiTraits();
	
	    /// destructor
	    virtual ~FeaFiTraits();
			
			/**
				@brief copy input data to external memory and update range information
				
				@p buffer_size is the size of the ring buffer used in the internal MSExperimentExtern
			*/
			template <class SpectrumIteratorType>
	    void setData(const SpectrumIteratorType& begin, const SpectrumIteratorType& end, UnsignedInt buffer_size)
		  {
		  	map_.setBufferSize( buffer_size );
				map_.updateBuffer();
		
				for (SpectrumIteratorType it = begin; it != end; ++it)
				{
					if (it->getMSLevel() == 1 && it->size() > 0) map_.push_back(*it);	// remove empty scans.
				}	
			
				std::cout << "Updating range information. " << std::endl;
		    // update range informations
		    map_.updateRanges();
		
				if (map_.getSize() == 0)
				{
					std::cout << "No data with MS level 1 provided. Aborting. " << std::endl;
					return;
				}
									
				std::cout << "This map contains " << map_.size() << " scans and " << map_.getSize() << " data points. " << std::endl;
		
		    // resize internal data structures
		    flags_.resize(map_.size());
				for (UnsignedInt i=0; i<map_.size(); ++i)
				{
					flags_[i].assign(map_[i].size(),FeaFiTraits::UNUSED);
				}
		  }

			
			/// Const access to LC-MS map
			inline const MapType& getData() const 
			{ 
				return map_; 
			}
				
	    /// non-mutable access flag with index @p index .
	    inline const Flag& getPeakFlag(const IDX& index) const
	    {
	    	return flags_[index.first][index.second];
	    }
	    /// mutable access flag with index @p index.
	    inline Flag& getPeakFlag(const IDX& index) 
	    { 
	    	return flags_[index.first][index.second];
	    }
		
	    /// access intensity of peak with index @p index.
	    inline const IntensityType& getPeakIntensity(const IDX& index) const
	    { 
	    	return map_[index.first][index.second].getIntensity(); 
	    }
	    /// access m/z of peak with index @p index .
	    inline const CoordinateType& getPeakMz(const IDX& index) const
	    { 
	    	return map_[index.first][index.second].getPos(); 
	    }
	    /// access retention time of peak with index @p index.
	    inline const CoordinateType& getPeakRt(const IDX& index) const
	    { 
	    	return map_[index.first].getRetentionTime();
	    }
	
	    /// returns the 2D coordinates of a peak (needed for models)
	    inline PositionType2D getPeakPos(const IDX& index) const
			{ 
				return PositionType2D(map_[index.first].getRetentionTime(),map_[index.first][index.second].getPos());
			}
	
	    /// fills @p index with the index of next peak in m/z dimension
	    inline void getNextMz(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
		  	//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
		    OPENMS_PRECONDITION(index.second<map_[index.first].size(), "Peak index outside of scan!");
    
	    	//At the last peak of this spectrum
	      if (index.second==map_[index.first].size()-1)
	      {
	      	throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index);
	      }
	
	      ++index.second;
	    }
	
	    /// fills @p index with the index of previous peak in m/z dimension
	    inline void getPrevMz(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
		  	//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
		    OPENMS_PRECONDITION(index.second<map_[index.first].size(), "Peak index outside of scan!");
    
	      //begin of scan
	      if (index.second==0)
	      {
	      	throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);
				}
				
	      --index.second;
	    }
	
	    /// fills @p index with the index of nearest peak in m/z dimension in the next spectrum
	    void getNextRt(IDX& index) throw (NoSuccessor, Exception::Precondition);
	
	    /// fills @p index with the index of nearest peak in m/z dimension in the previous spectrum
			void getPrevRt(IDX& index) throw (NoSuccessor, Exception::Precondition);
			
			//Calculates the convex hull of a index set and adds it to the feature
			void addConvexHull(const IndexSet& set, DFeature<2>& f) const;
	
	    /// run main loop
	    const DFeatureMap<2>& run(const std::vector<BaseSeeder*>& seeders,
	                              const std::vector<BaseExtender*>& extenders,
	                              const std::vector<BaseModelFitter*>& fitters);
	
		protected:
	  	/// Writes gnuplot output (only for debugging purposes)
	    void writeGnuPlotFile_(IndexSet peaks, bool last,int nr_feat);
	
	    /// Container for peak data
	    MapType map_;
	
	    /// Flags indicating whether a peak is unused, a seed or inside a feature region
	    std::vector< std::vector<Flag> > flags_;
	
	    /// The found features in the LC/MS map
	    DFeatureMap<2> features_;
	};

namespace Internal
{
	/// makes operator*() return intensity of the corresponding peak
	struct IntensityIterator 
		: FeaFiTraits::IndexSet::const_iterator
	{
    IntensityIterator ( FeaFiTraits::IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
      : FeaFiTraits::IndexSet::const_iterator(iter),
      	traits_(traits)
    {
    }
    
    FeaFiTraits::IntensityType operator * () const throw()
    {
    	return traits_->getPeakIntensity( FeaFiTraits::IndexSet::const_iterator::operator *() );
    }
    
		protected:
	    FeaFiTraits const * traits_;
	};
	
	/// Makes operator*() return mz of the corresponding peak
	struct MzIterator 
		: FeaFiTraits::IndexSet::const_iterator
	{
    MzIterator ( FeaFiTraits::IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
			: FeaFiTraits::IndexSet::const_iterator(iter), 
				traits_(traits)
    {
    }
    
    FeaFiTraits::CoordinateType operator * () const throw()
    {
    	return traits_->getPeakMz( FeaFiTraits::IndexSet::const_iterator::operator *() );
    }
    
		protected:
	    FeaFiTraits const * traits_;
	};
	
	/// Makes operator*() return retention time of the corresponding peak
	struct RtIterator 
		: FeaFiTraits::IndexSet::const_iterator
	{
    RtIterator ( FeaFiTraits::IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
    	: FeaFiTraits::IndexSet::const_iterator(iter),
    		traits_(traits)
    {
    }
    
    FeaFiTraits::CoordinateType operator * () const throw()
    {
    	return traits_->getPeakRt( FeaFiTraits::IndexSet::const_iterator::operator *() );
    }
		
		protected:
	    FeaFiTraits const * traits_;
};

} // namespace Internal

}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H
