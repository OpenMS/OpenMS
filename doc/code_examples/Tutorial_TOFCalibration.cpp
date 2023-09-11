// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;

  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);

  TOFCalibration ec;
  PeakMap exp_raw, calib_exp;
  MzMLFile mzml_file;
  mzml_file.load(tutorial_data_path + "/data/Tutorial_TOFCalibration_peak.mzML", calib_exp);
  mzml_file.load(tutorial_data_path + "/data/Tutorial_TOFCalibration_raw.mzML", exp_raw);

  vector<double> ref_masses;
  TextFile ref_file;
  ref_file.load(tutorial_data_path + "/data/Tutorial_TOFCalibration_masses.txt", true);
  for (TextFile::ConstIterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
  {
    ref_masses.push_back(String(iter->c_str()).toDouble());
  }

  std::vector<double> ml1;
  ml1.push_back(418327.924993827);

  std::vector<double> ml2;
  ml2.push_back(253.645187196031);

  std::vector<double> ml3;
  ml3.push_back(-0.0414243465397252);

  ec.setML1s(ml1);
  ec.setML2s(ml2);
  ec.setML3s(ml3);

  Param param;
  param.setValue("PeakPicker:peak_width", 0.1);
  ec.setParameters(param);
  ec.pickAndCalibrate(calib_exp, exp_raw, ref_masses);

  return 0;
} //end of main
