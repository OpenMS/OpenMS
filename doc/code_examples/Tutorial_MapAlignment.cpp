// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto tutorial_data_path = OPENMS_DOC_PATH + String("/code_examples/");

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
} // end of main
