// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [FeatureMap]

#include <OpenMS/KERNEL/FeatureMap.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // Insert of two features into a map and iterate over the features.
  FeatureMap map;

  Feature feature;
  feature.setRT(15.0);
  feature.setMZ(571.3);
  map.push_back(feature); //append feature 1
  feature.setRT(23.3);
  feature.setMZ(1311.3);
  map.push_back(feature); //append feature 2

  // Iteration over FeatureMap
  for (auto it = map.begin(); it != map.end(); ++it)
  {
    cout << it->getRT() << " - " << it->getMZ() << endl;
  }

  // Calculate and output the ranges
  map.updateRanges();
  cout << "Int: " << map.getMinIntensity() << " - " << map.getMaxIntensity() << endl;
  cout << "RT:  " << map.getMinRT() << " - " << map.getMaxRT() << endl;
  cout << "m/z: " << map.getMinMZ() << " - " << map.getMaxMZ() << endl;

  // ... and many more
  return 0;
} //end of main

//! [FeatureMap]
