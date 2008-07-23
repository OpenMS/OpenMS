#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap exp;
  
  MzDataFile mzdata_file;
  mzdata_file.load("data/Tutorial_GaussFilter.mzData",exp);

  GaussFilter g;
  Param param;
  param.setValue("gaussian_width",1.0);
  g.setParameters(param);
 
  g.filterExperiment(exp);
 
  return 0;
} //end of main
