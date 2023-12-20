// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

using namespace OpenMS;
using namespace std;

Int main()
{
  FeatureFinder ff;
  // ... set parameters (e.g. from INI file)
  Param parameters;
  // ... set input data (e.g. from mzML file)
  PeakMap input;
  // ... set output data structure
  FeatureMap output;
  // ... set user-specified seeds, if needed
  FeatureMap seeds;

  ff.run("simple", input, output, parameters, seeds);

  return 0;
} //end of main
