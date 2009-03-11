#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap exp;
  
  MzDataFile mzdata_file;
  mzdata_file.load("data/Tutorial_MorphologicalFilter.mzData",exp);

  Param parameters;
  parameters.setValue("struc_elem_length", 1.0);
  parameters.setValue("struc_elem_unit","Thomson");
  parameters.setValue("method","tophat");

  MorphologicalFilter th;
  th.setParameters(parameters);
 
  th.filterExperiment(exp);

  return 0;
} //end of main
