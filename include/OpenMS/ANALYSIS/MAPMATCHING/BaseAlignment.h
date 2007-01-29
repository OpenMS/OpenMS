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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEALIGNMENT_H

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringPairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairwiseMapMatcher_registerChildren.h>

#include <utility>
#include <fstream>
#include <vector>

namespace OpenMS
{
  /**
  Map alignment classes.

  This module contains all classes that can be used for the alignment of multiple maps.

  @defgroup MapAlignment Map Alignment
  */

  /**
     @brief Base alignment class.
     
     This class is the base class for the alignment of multiple element maps. 
     An element can be a DPeak, a DFeature, a ConsensusPeak or a ConsensusFeature.
     Corresponding elements are grouped together and stored as a consensus
     element in the final consensus map (a stl vector of consensus elements). 
   
  */
  template < typename ConsensusElementT >
  class BaseAlignment
  {
  public:
    /// Defines the coordinates of peaks / features.
    enum DimensionId
    {
      RT = DimensionDescription < LCMS_Tag >::RT,
      MZ = DimensionDescription < LCMS_Tag >::MZ
    };

    /// Consensus element type
    typedef ConsensusElementT ConsensusElementType;

    /// Element type (the type of the elements that are grouped together)
    typedef typename ConsensusElementType::ElementType ElementType;
    /// Element container type (the type of the element maps, that are aligned)
    typedef typename ConsensusElementType::ElementContainerType ElementContainerType;
    /// Consensus map type
    typedef std::vector< ConsensusElementType > ConsensusMapType;
    /// Type of estimated transformation
    typedef DGrid< 2 > GridType;

    /// Constructor
    BaseAlignment()
        : param_()
    {}

    /// Copy constructor 
    BaseAlignment(const BaseAlignment& source)
        : param_(source.param_),
        final_consensus_map_(source.final_consensus_map_),
        transformations_(source.transformations_),
        file_names_(source.file_names_),
        element_map_vector_(source.element_map_vector_),
        map_type_(source.map_type_)
    {}

    ///  Assignment operator 
    virtual BaseAlignment& operator = (const BaseAlignment& source)
    {
      if (&source==this)
        return *this;

      param_ = source.param_;
      final_consensus_map_ = source.final_consensus_map_;
      transformations_ = source.transformations_;
      file_names_ = source.file_names_;
      element_map_vector_ = source.element_map_vector_;
      map_type_ = source.map_type_;
      return *this;
    }

    /// Destructor
    virtual ~BaseAlignment()
  {}

    /// Get param class (non-mutable)
    const Param& getParam() const
    {
      return param_;
    }

    /// Mutable access to the param object
    void setParam(const Param& param)
    {
      param_ = param;
      map_type_ = param_.getValue("map_type");
    }

    /// Mutable access to the vector of maps
    void setElementMapVector(const std::vector< ElementContainerType* >& element_map_vector)
    {
      element_map_vector_ = element_map_vector;
    }
    /// Mutable access to the vector of maps
    std::vector< ElementContainerType* >& getElementMapVector()
    {
      return element_map_vector_;
    }
    /// Non-mutable access to the vector of maps
    const std::vector< ElementContainerType* >& getElementMapVector() const
    {
      return element_map_vector_;
    }

    /// Mutable access to the transformations
    void setTransformationVector(const std::vector< GridType >& transformations)
    {
      transformations_ = transformations;
    }
    /// Non-mutable access to the transformations
    const std::vector< GridType >& getTransformationVector() const
    {
      return transformations_;
    }

    /// Mutable access to the vector of file names
    void setFileNames(const std::vector< String >& file_names)
    {
      file_names_ = file_names;
    }
    /// Non-mutable access to the vector of file names
    const std::vector< String >& getFileNames() const
    {
      return file_names_;
    }

    /// Mutable access to the map type
    void setMapType(const String& map_type)
    {
      map_type_ = map_type;
      String param_name = "map_type";
      param_.setValue(param_name, map_type);
    }
    /// Non-mutable access to the map type (map type can be "feature_map", "peak_map" or "consensus_map"
    const String& getMapType() const
    {
      return map_type_;
    }

    /// Non-mutable access to the final consensus map
    const std::vector < ConsensusElementType >& getFinalConsensusMap() const
    {
      return final_consensus_map_;
    }
    /// Mutable access to the final consensus map
    void setFinalConsensusMap(const std::vector < ConsensusElementType >& final_consensus_map)
    {
      final_consensus_map_ = final_consensus_map;
    }

    /// Start the alignment
    virtual void run() throw (Exception::InvalidValue) = 0;


    /// Return the alignment tree in Newick Tree format
    virtual String getAlignmentTree() const = 0;

  protected:
    /// Parameter object
    Param param_;

    /// Final consensus map
    std::vector < ConsensusElementType > final_consensus_map_;

    /// The transformation vector
    std::vector< GridType > transformations_;

    /// File names of the maps to be aligned
    std::vector< String > file_names_;

    /// The maps to be aligned
    std::vector< ElementContainerType* > element_map_vector_;

    /// The map type
    String map_type_;

    /// Build a consensus map of the map with index map_index (the set of grouped elements contains the element itself).
    void buildConsensusMapTypeInsertInGroup_(UnsignedInt map_index, std::vector< ConsensusElementType >& cons_map)
    {
      const ElementContainerType& map = *(element_map_vector_[map_index]);
      UnsignedInt n = map.size();
      for (UnsignedInt i=0; i < n; ++i)
      {
        ConsensusElementType c(map_index,i,map[i]);
        cons_map.push_back(c);
      }
    }

    /// Build a consensus map of the map with index map_index (the set of grouped elements doesn't contain the element itself).
    void buildConsensusMapType_(UnsignedInt map_index, std::vector< ConsensusElementType >& cons_map)
    {
      const ElementContainerType& map = *(element_map_vector_[map_index]);
      UnsignedInt n = map.size();
      for (UnsignedInt i=0; i < n; ++i)
      {
        ConsensusElementType c(map[i].getPosition(),map[i].getIntensity());
        cons_map.push_back(c);
      }
    }
  }
  ; // BaseAlignment
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEALIGNMENT_H
