// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>

namespace OpenMS
{

  FeatureGroupingAlgorithmUnlabeled::FeatureGroupingAlgorithmUnlabeled() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmUnlabeled");
    defaults_.insert("", StablePairFinder().getParameters());
    defaultsToParam_();
    // The input for the pairfinder is a vector of FeatureMaps of size 2
    pairfinder_input_.resize(2);
  }

  FeatureGroupingAlgorithmUnlabeled::~FeatureGroupingAlgorithmUnlabeled() = default;

  void FeatureGroupingAlgorithmUnlabeled::group(const std::vector<FeatureMap> & maps, ConsensusMap & out)
  {
    // check that the number of maps is ok
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "At least two maps must be given!");
    }

    // define reference map (the one with most peaks)
    Size reference_map_index = 0;
    Size max_count = 0;
    for (Size m = 0; m < maps.size(); ++m)
    {
      if (maps[m].size() > max_count)
      {
        max_count = maps[m].size();
        reference_map_index = m;
      }
    }

    std::vector<ConsensusMap> input(2);

    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
    MapConversion::convert(reference_map_index, maps[reference_map_index],
                          input[0]);

    // loop over all other maps, extend the groups
    StablePairFinder pair_finder;
    pair_finder.setParameters(param_.copy("", true));

    for (Size i = 0; i < maps.size(); ++i)
    {
      if (i != reference_map_index)
      {
        MapConversion::convert(i, maps[i], input[1]);
        // compute the consensus of the reference map and map i
        ConsensusMap result;
        pair_finder.run(input, result);
        input[0].swap(result);
      }
    }

    // replace result with temporary map
    out.swap(input[0]);
    // copy back the input maps (they have been deleted while swapping)
    out.getColumnHeaders() = input[0].getColumnHeaders();

    postprocess_(maps, out);
  }

  void FeatureGroupingAlgorithmUnlabeled::addToGroup(int map_id, const FeatureMap& feature_map)
  {
    // create new PairFinder
    StablePairFinder pair_finder;
    pair_finder.setParameters(param_.copy("", true));

    // Convert the input map to a consensus map (using the given map_id) and
    // replace the second element in the pairfinder_input_ vector.
    MapConversion::convert(map_id, feature_map, pairfinder_input_[1]);

    // compute the consensus of the reference map and map map_id
    ConsensusMap result;
    pair_finder.run(pairfinder_input_, result);
    pairfinder_input_[0].swap(result);
  }

} // namespace OpenMS
