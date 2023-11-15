// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);

  PeakMap exp;

  FileHandler().loadExperiment(tutorial_data_path + "/data/Tutorial_GaussFilter.mzML", exp, {FileTypes::MZML});

  GaussFilter g;
  Param param;
  param.setValue("gaussian_width", 1.0);
  g.setParameters(param);

  g.filterExperiment(exp);

  return 0;
} //end of main
