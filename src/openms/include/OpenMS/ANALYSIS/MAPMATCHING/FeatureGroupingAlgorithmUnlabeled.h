// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

#include <OpenMS/KERNEL/ConversionHelper.h>

namespace OpenMS
{
  /**
      @brief A map feature grouping algorithm for unlabeled data.

      It takes multiple maps and searches for corresponding features.
      The corresponding features must be aligned, but may have small position deviations.

      There are two ways to run the algorithm:

      - a) Simply call "group" with all maps in memory.
      - b) Call "setReference", "addToGroup" (n times), "getResultMap" in that order.
      
      The second way is more memory efficient because at all times, only the
      reference map and the current map need to be in memory

      @htmlinclude OpenMS_FeatureGroupingAlgorithmUnlabeled.parameters

      @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI FeatureGroupingAlgorithmUnlabeled :
    public FeatureGroupingAlgorithm
  {
public:
    /// Default constructor
    FeatureGroupingAlgorithmUnlabeled();

    /// Destructor
    ~FeatureGroupingAlgorithmUnlabeled() override;

    /**
        @brief Sets the reference map for the algorithm
        
        This will store the input map in the first argument of the pairfinder_input_ vector.

        Note that map_id needs to be unique (see addToGroup docu).
    */
    template <typename MapType>
    void setReference(int map_id, const MapType & map)
    {
      MapConversion::convert(map_id, map, pairfinder_input_[0]);
    }

    /**
        @brief Returns the computed consensus map (after calling addToGroup with all maps)
    */
    ConsensusMap & getResultMap()
    {
      return pairfinder_input_[0];
    }

    /**
        @brief Applies the algorithm

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    void group(const std::vector<FeatureMap > & maps, ConsensusMap & out) override;

    /**
        @brief Adds one map to the group

        This assumes that setReference has been called before. It is the callers
        responsibility to provide a valid map_id. Note that all calls to
        addToGroup _must_ use different map_ids and they all need to be
        different from the one used to call setReference!

        @exception IllegalArgument is thrown if the map_id has been used before
    */
    void addToGroup(int map_id, const FeatureMap & feature_map);

private:

    // This vector should always have 2 elements
    // - the first element is the currently computed consensus map.
    //   After initialization of the algorithm, it will consist of the reference
    //   map alone, after adding all maps through addToGroup it will contain the
    //   final result (e.g. the consensus result of all maps).
    // - the second element is the map which was last added to the consensus map
    std::vector<ConsensusMap> pairfinder_input_;

    /// Copy constructor intentionally not implemented -> private
    FeatureGroupingAlgorithmUnlabeled(const FeatureGroupingAlgorithmUnlabeled &);
    /// Assignment operator intentionally not implemented -> private
    FeatureGroupingAlgorithmUnlabeled & operator=(const FeatureGroupingAlgorithmUnlabeled &);

  };

} // namespace OpenMS

