// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  FeatureMap reference;
  FeatureMap toAlign;

  FileHandler xml_file;
  xml_file.loadFeatures(tutorial_data_path + "/data/Tutorial_MapAlignment_1.featureXML", reference);
  xml_file.loadFeatures(tutorial_data_path + "/data/Tutorial_MapAlignment_2.featureXML", toAlign);

  // create map alignment algorithm
  MapAlignmentAlgorithmPoseClustering algorithm;
  
  // ... set parameters
  algorithm.setReference(reference);

  // create object for the computed transformation
  TransformationDescription transformation;

  // align
  algorithm.align(toAlign, transformation);

  // store results
  xml_file.storeFeatures("Tutorial_MapAlignment_1.featureXML", reference);
  xml_file.storeFeatures("Tutorial_MapAlignment_2.featureXML", toAlign);

  return 0;
} //end of main
