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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
    @brief The base class of all superimposer algorithms.

     This class defines the basic interface for all superimposer
     algorithms. It works on two element maps and
     computes a transformation, that maps the elements of one map (scene map) 
     as near as possible to the elements in the other map (model map).
     
     The element map must be a random access container (e.g. vector, DPeakArray, FeatureMap)
     of elements that have the same interface as RawDataPoint2D.
  */
  template <typename MapT>
  class BaseSuperimposer 
  	: public FactoryProduct
  {
  	
  	public:

	    /// Container for input elements
	    typedef MapT ElementMapType;
	
	    /// Constructor
	    BaseSuperimposer()
	    	: FactoryProduct("BaseSuperimposer"),
          scene_map_(0),
          model_map_(0)
	    {
	    }
	
	    /// Destructor
	    virtual ~BaseSuperimposer()
	  	{
	  	}
	  	
	    /// Set  model map (stores a pointer to the map only, so it must still exist, when calling run())
	    virtual void setModelMap(const ElementMapType& map)
	    {
	      model_map_ = &map;
	    }
	
	    /// Sets the scene map (stores a pointer to the map only, so it must still exist, when calling run())
	    virtual void setSceneMap(const ElementMapType& map)
	    {
	      scene_map_ = &map;
	    }
	
	    /**
	    	@brief Estimates the transformation and fills the given mapping function
	    	
	    	@exception IllegalArgument is thrown if the element maps are not set.
	    */
	    virtual void run(LinearMapping& mapping) = 0;
	
	    /// Register all derived classes here
	    static void registerChildren();
	
	  protected:
	    /// scene map
	    const ElementMapType* scene_map_;
	    /// model map
			const ElementMapType* model_map_;
  
  }; // BaseSuperimposer

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H
