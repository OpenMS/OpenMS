#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  
  MzXMLFile mzxml;
  MzMLFile mzml;

  // temporary data storage
  MSExperiment<Peak1D> map;

  // convert MzXML to MzML
  mzxml.load("data/Tutorial_FileIO.mzXML",map);
  mzml.store("output/Tutorial_FileIO.mzML",map);
  
  return 0;
} //end of main
