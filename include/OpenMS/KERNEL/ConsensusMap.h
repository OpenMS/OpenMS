// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/KERNEL/RangeManager.h>

namespace OpenMS
{
	/**
    @brief A container for consensus elements.
    
    A %ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as a vector of elements.
    
    The map indices used in the consensus features should be registered in this class.
 		
 		@improvement Add method to dump gnuplot files of consensus maps (Clemens)
 		
    @ingroup Kernel
  */
	class ConsensusMap 
		: public DPeakArray<ConsensusFeature>,
			public RangeManager<2>
	{
	  public:
	  	/// Source file desciption for input files
	  	struct FileDescription
	  		: public MetaInfoInterface
	  	{
	  		///Default constructor
	  		FileDescription()
	  			: MetaInfoInterface(),
	  				filename(),
	  				label(),
	  				size(0)
	  		{
	  		}

	  		/// file name of the file
	  		String filename;
	  		/// Label e.g. 'heavy' and 'light' for ICAT, or 'sample1' and 'sample2' for label-free quantitation
	  		String label;
	  		/// @brief Number of elements (features, peaks, ...).
	  		/// This is e.g. used to check for correct element indices when writing a consensus map
	  		UInt size;
	  	};

	    ///@name Type definitions
	    //@{
	    typedef DPeakArray<ConsensusFeature > Base;
			typedef RangeManager<2> RangeManagerType;
	  	typedef Map<UInt,FileDescription> FileDescriptions;
			//@}
	  	
	    /// Default onstructor
	    inline ConsensusMap()
	      : Base(),
	      	RangeManagerType()
	    {
	    }
	
	    /// Copy constructor
	    inline ConsensusMap(const ConsensusMap& source)
	      : Base(source),
	      	RangeManagerType(source),
	        file_description_(source.file_description_)
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
	    ConsensusMap& operator=(const ConsensusMap& source)
	    {
	      if (this==&source) return *this;
	
	      Base::operator=(source);
				RangeManagerType::operator=(source);
	      file_description_ = source.file_description_;
	      
	      return *this;
	    }
	
	    /// Non-mutable access to the file descriptions
	    inline const FileDescriptions& getFileDescriptions() const
	    {
	      return file_description_;
	    }

	    /// Mutable access to the file descriptions
	    inline FileDescriptions& getFileDescriptions()
	    {
	      return file_description_;
	    }
	    			
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
				
				output_map.getFileDescriptions()[input_map_index].size = input_map.size();
			}
			
			/**
				@brief Similar to convert, but copies only the @p n most intense elements from an MSExperiment.
		
				@param input_map_index The index of the input map.
				@param input_map The input map to be converted.
				@param output_map The resulting ConsensusMap.
				@param n The maximum number of elements to be copied.
			*/
			static void convert( UInt const input_map_index, MSExperiment<> & input_map, ConsensusMap& output_map, UInt n )
			{
				input_map.updateRanges(1);
				if ( n > input_map.getSize() )
				{
					n = input_map.getSize();
				}
				output_map.clear();
				output_map.reserve(n);
				std::vector<Peak2D> tmp;
				tmp.reserve(input_map.getSize());
				input_map.get2DData(tmp); //Avoid tripling the memory consumption by this call
				std::partial_sort( tmp.begin(), tmp.begin()+n, tmp.end(), reverseComparator(Peak2D::IntensityLess()) );
				for ( UInt element_index = 0; element_index < n; ++element_index )
				{
					output_map.push_back( ConsensusFeature( input_map_index, element_index, tmp[element_index] ) );
				}
				
				output_map.getFileDescriptions()[input_map_index].size = n;
			}
			
			// Docu in base class
			void updateRanges();
			
	  protected:
	    /// Map from index to file description
	  	FileDescriptions file_description_;
  };

  ///Print the contents of a ConsensusMap to a stream.
  std::ostream& operator << (std::ostream& os, const ConsensusMap& cons_map);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSMAP_H
