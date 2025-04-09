// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/KERNEL/FeatureMap.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  FeatureMap map;

  Feature feature;
  feature.setIntensity(461.3f);
  feature.setRT(15.0);
  feature.setMZ(571.3);
  map.push_back(feature);
  feature.setIntensity(12213.5f);
  feature.setRT(23.3);
  feature.setMZ(1311.3);
  map.push_back(feature);

  //calculate the ranges
  map.updateRanges();

  cout << "Int: " << map.getMinIntensity() << " - " << map.getMaxIntensity() << endl;
  cout << "RT:  " << map.getMinRT() << " - " << map.getMaxRT() << endl;
  cout << "m/z: " << map.getMinMZ() << " - " << map.getMaxMZ() << endl;

  return 0;
} //end of main
