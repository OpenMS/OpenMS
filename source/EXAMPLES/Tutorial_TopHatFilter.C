#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap exp_raw;
  PeakMap exp_filtered;
  
  MzDataFile mzdata_file;
  mzdata_file.load("data/Tutorial_TopHatFilter.mzData",exp_raw);

  TopHatFilter th;
  Param param;
  param.setValue("struc_elem_length",1.0);
  th.setParameters(param);
 
  th.filterExperiment(exp_raw,exp_filtered);
 
  return 0;
} //end of main
