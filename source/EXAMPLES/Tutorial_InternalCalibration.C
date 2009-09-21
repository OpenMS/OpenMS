#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  InternalCalibration ic;
  PeakMap exp_raw,exp_calibrated;
  MzMLFile mzml_file;
  mzml_file.load("data/Tutorial_InternalCalibration.mzML",exp_raw);

  std::vector<double> ref_masses;
  ref_masses.push_back(1296.68476942);
  ref_masses.push_back(2465.19833942);

  Param param;
  param.setValue("PeakPicker:thresholds:peak_bound",800.0);
  param.setValue("PeakPicker:peak_width",0.15);
  ic.setParameters(param);
  
  ic.calibrateMapSpectrumwise(exp_raw,exp_calibrated,ref_masses);
  
  return 0;
} //end of main
