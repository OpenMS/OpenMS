// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto tutorial_data_path = OPENMS_DOC_PATH + String("/code_examples/");
  
  vector<FeatureMap > maps;
  maps.resize(2);

  FileHandler feature_file;
  feature_file.loadFeatures(tutorial_data_path + "/data/Tutorial_Unlabeled_1.featureXML", maps[0]);
  feature_file.loadFeatures(tutorial_data_path + "/data/Tutorial_Unlabeled_2.featureXML", maps[1]);

  ConsensusMap out;
  out.getColumnHeaders()[0].filename = "/data/Tutorial_Unlabeled_1.mzML";
  out.getColumnHeaders()[0].size = maps[0].size();
  out.getColumnHeaders()[1].filename = "/data/Tutorial_Unlabeled_2.mzML";
  out.getColumnHeaders()[1].size = maps[1].size();


  FeatureGroupingAlgorithmUnlabeled algorithm;
  // ... set parameters
  algorithm.group(maps, out);
  FileHandler consensus_file;
  consensus_file.storeConsensusFeatures("Tutorial_Unlabeled.consensusXML", out);

  return 0;
} //end of main
