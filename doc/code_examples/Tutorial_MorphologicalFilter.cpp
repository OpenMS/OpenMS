// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/openms_data_path.h> // exotic header for path to tutorial data
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  auto tutorial_data_path = OPENMS_DOC_PATH + String("/code_examples/");

  PeakMap exp;

  FileHandler().loadExperiment(tutorial_data_path + "/data/Tutorial_MorphologicalFilter.mzML", exp);

  Param parameters;
  parameters.setValue("struc_elem_length", 1.0);
  parameters.setValue("struc_elem_unit", "Thomson");
  parameters.setValue("method", "tophat");

  MorphologicalFilter mf;
  mf.setParameters(parameters);

  mf.filterExperiment(exp);

  return 0;
} // end of main
