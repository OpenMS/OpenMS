#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  RawSpectrum spec_raw;
  RawSpectrum spec_filtered;
  
  DTAFile dta_file;
  dta_file.load("../TEST/data/PeakTypeEstimator_raw.dta",spec_raw);

  SavitzkyGolaySVDFilter sg;
  sg.setWindowSize(21);
  sg.setOrder(3);
  
  sg.filter(spec_raw,spec_filtered);
 
  return 0;
} //end of main
