#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  RawMap exp_raw;
  PeakMap exp_picked;
  
  MzDataFile mzdata_file;
  mzdata_file.load("../TEST/data/PeakPicker_test.mzData",exp_raw);

  PeakPickerCWT pp;
  pp.setWaveletScale(0.2);
  pp.setFwhmBound(0.1);
  pp.setPeakBound(500);

  pp.pickExperiment(exp_raw,exp_picked);
  exp_picked.updateRanges();
  
  cout << "Scale of the wavelet: " << pp.getWaveletScale()
  		 << "\nMinimal fwhm of a mass spectrometric peak: " << pp.getFwhmBound()
  		 << "\nMinimal intensity of a mass spectrometric peak " << pp.getPeakBound()
  		 << "\n\nNumber of picked peaks " << exp_picked.getSize() << std::endl;

  return 0;
} //end of main
