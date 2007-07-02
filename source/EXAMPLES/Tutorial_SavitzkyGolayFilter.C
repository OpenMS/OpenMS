#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  RawSpectrum spec_raw;
  RawSpectrum spec_resampled;
  RawSpectrum spec_filtered;
  
  DTAFile dta_file;
  dta_file.load("../TEST/data/PeakTypeEstimator_rawTOF.dta",spec_raw);
  
  LinearResampler lr;
  lr.setSpacing(0.01);
  lr.raster(spec_raw,spec_resampled);

  SavitzkyGolaySVDFilter sg;
  sg.setWindowSize(21);
  sg.setOrder(3);
  
  sg.filter(spec_resampled,spec_filtered);
 
  return 0;
} //end of main
