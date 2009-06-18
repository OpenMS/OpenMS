#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap exp;
  
  MzMLFile mzml_file;
  mzml_file.load("data/Tutorial_MorphologicalFilter.mzML",exp);

  Param parameters;
  parameters.setValue("struc_elem_length", 1.0);
  parameters.setValue("struc_elem_unit","Thomson");
  parameters.setValue("method","tophat");

  MorphologicalFilter th;
  th.setParameters(parameters);
 
  th.filterExperiment(exp);

  return 0;
} //end of main
