// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  vector<FeatureMap> maps(1, FeatureMap{});

  FileHandler().loadFeatures(OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_Labeled.featureXML"), maps[0], {FileTypes::FEATUREXML});
  ConsensusMap out;
  out.getColumnHeaders()[0].filename = "data/Tutorial_Labeled.mzML";
  out.getColumnHeaders()[0].size = maps[0].size();
  out.getColumnHeaders()[0].label = "light";
  out.getColumnHeaders()[1].filename = "data/Tutorial_Labeled.mzML";
  out.getColumnHeaders()[1].size = maps[0].size();
  out.getColumnHeaders()[1].label = "heavy";

  FeatureGroupingAlgorithmLabeled algorithm;
  // ... set parameters
  algorithm.group(maps, out);
  FileHandler().storeConsensusFeatures("Tutorial_Labeled.consensusXML", out, {FileTypes::CONSENSUSXML});

  return 0;
} //end of main
