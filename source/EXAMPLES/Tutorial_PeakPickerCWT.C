#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap exp_raw;
  PeakMap exp_picked;
  
  MzDataFile mzdata_file;
  mzdata_file.load("data/Tutorial_PeakPickerCWT.mzData",exp_raw);

  PeakPickerCWT pp;
  Param param;
  param.setValue("thresholds:peak_bound",500.0);
  param.setValue("thresholds:fwhm_bound",0.1);
  param.setValue("wavelet_transform:scale",0.2);
  pp.setParameters(param);

  pp.pickExperiment(exp_raw,exp_picked);
  exp_picked.updateRanges();
  
  cout << "Scale of the wavelet: " << (DoubleReal)param.getValue("wavelet_transform:scale")
       << "\nMinimal fwhm of a mass spectrometric peak: " << (DoubleReal)param.getValue("thresholds:fwhm_bound")
       << "\nMinimal intensity of a mass spectrometric peak " << (DoubleReal)param.getValue("thresholds:peak_bound")
       << "\n\nNumber of picked peaks " << exp_picked.getSize() << std::endl;

  return 0;
} //end of main
