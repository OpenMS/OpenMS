#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  TOFCalibration ec;
  PeakMap exp_raw,calib_exp;
  MzMLFile mzml_file;
  mzml_file.load("data/Tutorial_TOFCalibration_peak.mzML",calib_exp);
  mzml_file.load("data/Tutorial_TOFCalibration_raw.mzML",exp_raw);
  
  vector<DoubleReal> ref_masses;
  TextFile ref_file;
  ref_file.load("data/Tutorial_TOFCalibration_masses.txt",true);
  for(TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
  {
    ref_masses.push_back(String(iter->c_str()).toDouble());
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
  param.setValue("PeakPicker:peak_width",0.1);
  ec.setParameters(param);
  ec.pickAndCalibrate(calib_exp,exp_raw,ref_masses);

  return 0;
} //end of main
