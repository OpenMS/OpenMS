#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap exp_raw;
  PeakMap exp_picked;
  
  MzMLFile mzml_file;
  mzml_file.load("data/Tutorial_PeakPickerCWT.mzML",exp_raw);

  PeakPickerCWT pp;
  Param param;
  param.setValue("peak_width",0.1);
  pp.setParameters(param);

  pp.pickExperiment(exp_raw,exp_picked);
  exp_picked.updateRanges();
  
  cout << "\nMinimal fwhm of a mass spectrometric peak: " << (DoubleReal)param.getValue("peak_width")
       << "\n\nNumber of picked peaks " << exp_picked.getSize() << std::endl;

  return 0;
} //end of main
