#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  
  MzXMLFile mzxml;
  MzDataFile mzdata;

  // temporary data storage
  MSExperiment<Peak1D> map;

  // convert MzXML to MzData
  mzxml.load("data/Tutorial_FileIO.mzXML",map);
  mzdata.store("output/Tutorial_FileIO.mzData",map);
  
  return 0;
} //end of main
