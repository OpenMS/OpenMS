// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{
	/**
    @brief A container for consensus elements.
    
    A ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as a vector of elements.
    
    The map indices used in the consensus features should be registered in this class.
 		
 		@todo Add method to dump gnuplot files of consensus maps (Clemens)
 		
    @ingroup Kernel
  */
	class ConsensusMap : public DPeakArray<ConsensusFeature>
	{
	  public:
	    /// Base class type
	    typedef DPeakArray<ConsensusFeature > Base;
	  	
	    /// Default onstructor
	    inline ConsensusMap()
	      : Base()
	    {
	    }
	
	    /// Copy constructor
	    inline ConsensusMap(const ConsensusMap& source)
	      : Base(source),
	        filenames_(source.filenames_)
	    {
	    }
	
	    /// Destructor
	    inline ~ConsensusMap()
	    {
	    }
	
	    /// Creates a ConsensusMap with n elements
	    inline ConsensusMap(Base::size_type n) 
	    	: Base(n)
	    {
	    }
	
	    /// Assignment operator
	    ConsensusMap& operator = (const ConsensusMap& source)
	    {
	      if (this==&source) return *this;
	
	      Base::operator=(source);
	      filenames_ = source.filenames_;
	      
	      return *this;
	    }
	
	    /// Non-mutable access to the filenames
	    inline const Map<UInt,String>& getFileNames() const
	    {
	      return filenames_;
	    }
	    
	    /// Set a file name
	    inline void setFileName(UInt index, const String& name)
	    {
	      filenames_[index] = name;
	    }
	    
	    /// Merge overlapping consensus elements
	    void merge(ConsensusMap& new_map);
			
			///Checks if all map identifiers in FeatureHandles are have a filename associated
			bool isValid() const;
			
			/// Sorts the peaks according to ascending quality
			void sortByQuality(bool reverse=false) 
			{ 
				if (reverse)
				{
					std::sort(Base::begin(), Base::end(), reverseComparator(ConsensusFeature::QualityLess())); 
				}
				else
				{
					std::sort(Base::begin(), Base::end(), ConsensusFeature::QualityLess()); 
				}
			}

			/**
				@brief Convert any (random access) container of features to a ConsensusMap.  Each
				ConsensusFeature contains a map index, so this has to be given as well.
				The previous content of output_map is cleared.
		
				@param input_map_index The index of the input map.
				@param input_map The container to be converted.  (Must support size() and operator[].)
				@param output_map The resulting ConsensusMap.
			*/
			template <typename ContainerT>
			static void convert( UInt const input_map_index, ContainerT const & input_map, ConsensusMap& output_map )
			{
				output_map.clear();
				output_map.reserve(input_map.size());
				for ( UInt element_index = 0; element_index < input_map.size(); ++element_index )
				{
					output_map.push_back( ConsensusFeature( input_map_index, element_index, input_map[element_index] ) );
				}
				return;
			}
			
			/**
				@brief Similar to convert, but copies only the @p n most intense elements from an MSExperiment.
		
				@param input_map_index The index of the input map.
				@param input_map The input map to be converted.
				@param output_map The resulting ConsensusMap.
				@param n The maximum number of elements to be copied.
			*/
			static void convert( UInt const input_map_index, MSExperiment<> & input_map, ConsensusMap& output_map, UInt n ) // TODO find out what goes wrong in template instantiation (?!!)
			// template <typename PeakT, typename AllocT>
			// static void convert( UInt const input_map_index, MSExperiment<PeakT,AllocT> & input_map, ConsensusMap& output_map, UInt n )
			{
				input_map.updateRanges(1);
				if ( n > input_map.getSize() )
				{
					n = input_map.getSize();
				}
				output_map.clear();
				output_map.reserve(n);
				std::vector<RawDataPoint2D> tmp; // TODO let's see if this will pass the nightly build
				// std::vector<RawDataPoint2D,AllocT> tmp;
				tmp.reserve(input_map.getSize());
				input_map.get2DData(tmp);
				std::partial_sort( tmp.begin(), tmp.begin()+n, tmp.end(), reverseComparator(RawDataPoint2D::IntensityLess()) );
				for ( UInt element_index = 0; element_index < n; ++element_index )
				{
					output_map.push_back( ConsensusFeature( input_map_index, element_index, tmp[element_index] ) );
				}
				return;
			}
			
			
	  protected:
	    /// Map from index to filenames
	    /// @todo Make filenames_ a map<UInt, STRUCT> with STRUCT = MetaInfoInterface, filename, what else? (Marc, Clemens)
			Map<UInt,String> filenames_;
  };

  ///Print the contents of a ConsensusMap to a stream.
  std::ostream& operator << (std::ostream& os, const ConsensusMap& cons_map);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSMAP_H
