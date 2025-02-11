// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>

#include <OpenMS/KERNEL/ConversionHelper.h>

namespace OpenMS
{

  FeatureGroupingAlgorithmLabeled::FeatureGroupingAlgorithmLabeled() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmLabeled");

    defaults_.insert("", LabeledPairFinder().getParameters());

    defaultsToParam_();
  }

  FeatureGroupingAlgorithmLabeled::~FeatureGroupingAlgorithmLabeled() = default;

  void FeatureGroupingAlgorithmLabeled::group(const std::vector<FeatureMap> & maps, ConsensusMap & out)
  {
    //check that the number of maps is ok
    if (maps.size() != 1)
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Exactly one map must be given!");
    if (out.getColumnHeaders().size() != 2)
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Two file descriptions must be set in 'out'!");

    //initialize LabeledPairFinder
    LabeledPairFinder pm;
    pm.setParameters(param_.copy("", true));

    //convert to consensus map
    std::vector<ConsensusMap> input(1);
    MapConversion::convert(0, maps[0], input[0]);

    //run
    pm.run(input, out);
  }

} //namespace
