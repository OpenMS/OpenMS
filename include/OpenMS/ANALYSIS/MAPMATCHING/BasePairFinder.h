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
		algorithms. It works on two consensus maps.
		
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
			virtual void setModelMap(ConsensusMap const& model_map)
			{
				model_map_ = &model_map;
			}
	
			/// Get model map
			ConsensusMap const & getModelMap() const
			{
				return *model_map_;
			}
	
			/// Set scene map.  @sa setModelMap()
			virtual void setSceneMap(ConsensusMap const& scene_map)
			{
				scene_map_ = &scene_map;
			}
	
			/// Get scene map
			ConsensusMap const & getSceneMap() const
			{
				return *scene_map_;
			}
	
			/// Run the algorithm
			virtual void run(ConsensusMap&)
			{
			}
	
	    /// Register all derived classes here
	    static void registerChildren();
	
	  protected:
			///@name map pointers and indices
			//@{
			const ConsensusMap* model_map_;
			const ConsensusMap* scene_map_;
			//@}
			
		 private:
	
	    /// Copy constructor intentionally not implemented
	    BasePairFinder(const BasePairFinder&);
			
	    /// Assignment operator intentionally not implemented
	    BasePairFinder & operator=(const BasePairFinder&);
		
	};

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
