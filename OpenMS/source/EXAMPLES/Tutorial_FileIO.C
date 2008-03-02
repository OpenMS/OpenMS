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
	MSExperiment<RawDataPoint1D> map;

	// convert MzXML to MzData
	mzxml.load("Tutorial_FileIO.mzXML",map);
	mzdata.store("Tutorial_FileIO.mzData",map);
	
  return 0;
} //end of main
