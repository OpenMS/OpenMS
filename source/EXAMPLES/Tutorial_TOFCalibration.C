#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  TOFCalibration ec;
  PeakMap exp_raw,calib_exp;
  MzDataFile mzdata_file;
  mzdata_file.load("data/Tutorial_TOFCalibration_peak.mzData",calib_exp);
  mzdata_file.load("data/Tutorial_TOFCalibration_raw.mzData",exp_raw);
  
  vector<DoubleReal> ref_masses;
  TextFile ref_file;
  ref_file.load("data/Tutorial_TOFCalibration_masses.txt",true);
  for(TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
  {
    ref_masses.push_back(atof(iter->c_str()));
  }

  std::vector<DoubleReal> ml1;
  ml1.push_back(418327.924993827);

  std::vector<DoubleReal> ml2;
  ml2.push_back(253.645187196031);
 
  std::vector<DoubleReal> ml3;
  ml3.push_back(-0.0414243465397252);

  ec.setML1s(ml1);
  ec.setML2s(ml2);
  ec.setML3s(ml3);

  Param param;
  param.setValue("PeakPicker:thresholds:peak_bound",800.0);
  param.setValue("PeakPicker:thresholds:fwhm_bound",0.1);
  param.setValue("PeakPicker:wavelet_transform:scale",0.12);
  ec.setParameters(param);
  ec.pickAndCalibrate(calib_exp,exp_raw,ref_masses);

  return 0;
} //end of main
