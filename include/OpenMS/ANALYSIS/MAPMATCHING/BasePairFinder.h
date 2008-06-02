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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{
 
  /**
		@brief The base class of all element pair finding algorithms.
			
		This class defines the basic interface for all element pair finding
		algorithms. It works on two element maps (FeatureMap is the default map
		type) and a transformation defined for the second element map (if no
		transformation is given, the pairs are found in the two original maps).  A
		element can be a DPeak, a DFeature or ConsensusFeature (wheras DFeature is
		the default element type).
  */
  class BasePairFinder
  	: public FactoryProduct
  {
	  public:
	    /// Default constructor
	    BasePairFinder();
	
			/// Destructor
	    virtual ~BasePairFinder();
	
			/** 
				@brief Set model map (might do some preprocessing).
	
				@param map_index If map_index>=0, then run() will use it as a map index
				for the input and store feature handles pointing to the consensus features
				of theq input to the consensus features in the result, which adds another
				level of indirection/nesting.<br> If -1, then run() will "unpack" the
				consensus features from the input, i.e. it will store the feature handles
				contained in the consensus features rather than the consensus features
				themselves in the result.
		
				@param model_map Consensus map to be used as model.
			*/
			virtual void setModelMap(Int map_index, ConsensusMap const& model_map)
			{
				model_map_ = &model_map;
				model_index_ = map_index;
			}
	
			/// Get model map
			ConsensusMap const & getModelMap() const
			{
				return *model_map_;
			}
	
			/// Set scene map.  @sa setModelMap()
			virtual void setSceneMap(Int map_index, ConsensusMap const& scene_map)
			{
				scene_map_ = &scene_map;
				scene_index_ = map_index;
			}
	
			/// Get scene map
			ConsensusMap const & getSceneMap() const
			{
				return *scene_map_;
			}
	
			/// Run the algorithm
			virtual void run(ConsensusMap& /*result_map*/)
			{
			};
	
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
	
	    /// Register all derived classes here
	    static void registerChildren();
	
	  protected:
			///@name map pointers and indices
			//@{
			const ConsensusMap* model_map_;
			const ConsensusMap* scene_map_;
			ConsensusMap* result_map_;
			Int model_index_;
			Int scene_index_;
			//@}
			
		 private:
	
	    /// Copy constructor intentionally not implemented
	    BasePairFinder(const BasePairFinder&);
			
	    /// Assignment operator intentionally not implemented
	    BasePairFinder & operator=(const BasePairFinder&);
		
	};

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
