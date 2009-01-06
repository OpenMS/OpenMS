#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap ms_exp;
  
  MzDataFile mzdata_file;
  mzdata_file.load("data/Tutorial_MorphologicalFilter.mzData",ms_exp);

  MorphologicalFilter th;
 
  th.filterMSExperiment(MorphologicalFilter::TOPHAT, 1.0, true, ms_exp);

  return 0;
} //end of main
