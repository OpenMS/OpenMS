// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  vector<FeatureMap > maps;
  maps.resize(1);

  FileHandler().loadFeatures(tutorial_data_path + "/data/Tutorial_Labeled.featureXML", maps[0], {FileTypes::FEATUREXML});
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
