// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>

using namespace OpenMS;
using namespace std;

Int main()
{
  FeatureFinderAlgorithmPicked ff;
  // ... set parameters (e.g. from INI file)
  Param parameters;
  // ... set input data (e.g. from mzML file)
  PeakMap input;
  // ... set output data structure
  FeatureMap output;
  // ... set user-specified seeds, if needed
  FeatureMap seeds;

  ff.run(input, output, parameters, seeds);

} //end of main
