// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmWNet.h>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithmWNet::FeatureGroupingAlgorithmWNet() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmWNet");

    defaults_.setValue("distance_metric", "LINF", "Distance metric for comparing feature positions");
    defaults_.setValidStrings("distance_metric", {"L1", "L2", "LINF"});

    defaults_.setValue("max_distance", 800.0, "Maximum distance between features to consider a match");
    defaults_.setMinFloat("max_distance", 0.0);

    defaults_.setValue("trash_cost", 800.0, "Cost of leaving a feature unmatched (higher = prefer matching over leaving unmatched)");
    defaults_.setMinFloat("trash_cost", 0.0);

    defaults_.setValue("mz_scale_factor", 0.0,
      "Scale factor applied to m/z values before distance computation. "
      "If 0 (default), auto-computed as max_RT_range / max_MZ_range to bring "
      "both dimensions to comparable numeric scale.");
    defaults_.setMinFloat("mz_scale_factor", 0.0);

    defaults_.setValue("normalize_intensities", "true",
      "Normalize feature intensities per map before alignment");
    defaults_.setValidStrings("normalize_intensities", {"true", "false"});

    defaultsToParam_();
  }

  FeatureGroupingAlgorithmWNet::~FeatureGroupingAlgorithmWNet() = default;

  template <typename MapType>
  void FeatureGroupingAlgorithmWNet::group_(const vector<MapType>& maps,
                                            ConsensusMap& out)
  {
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "At least two maps must be given!");
    }

    // TODO: implement pairwise WNetAligner<2> alignment and consensus building

    postprocess_(maps, out);
  }

  void FeatureGroupingAlgorithmWNet::group(const vector<FeatureMap>& maps,
                                           ConsensusMap& out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmWNet::group(const vector<ConsensusMap>& maps,
                                           ConsensusMap& out)
  {
    group_(maps, out);
  }

} // namespace OpenMS
