// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto file_gauss = OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_GaussFilter.mzML");

  PeakMap exp;

  FileHandler().loadExperiment(file_gauss, exp, {FileTypes::MZML});

  GaussFilter g;
  Param param;
  param.setValue("gaussian_width", 1.0);
  g.setParameters(param);

  g.filterExperiment(exp);

  return 0;
} //end of main
