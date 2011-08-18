#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  InternalCalibration ic;
  PeakMap exp,exp_calibrated;
  MzMLFile mzml_file;
  mzml_file.load("data/Tutorial_InternalCalibration.mzML",exp);

  std::vector<double> ref_masses;
  ref_masses.push_back(1296.68476942);
  ref_masses.push_back(2465.19833942);

  ic.calibrateMapSpectrumwise(exp,exp_calibrated,ref_masses);
  
  return 0;
} //end of main
