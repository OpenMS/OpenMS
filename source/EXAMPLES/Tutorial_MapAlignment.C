#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

using namespace OpenMS;
using namespace std;

Int main()
{
	vector<FeatureMap<> > maps;
	maps.resize(2);

  FeatureXMLFile xml_file;
  xml_file.load("data/Tutorial_MapAlignment_1.featureXML",maps[0]);
  xml_file.load("data/Tutorial_MapAlignment_2.featureXML",maps[1]);
  
  MapAlignmentAlgorithmPoseClustering algorithm;
  algorithm.alignFeatureMaps(maps);
  
  xml_file.store("output/Tutorial_MapAlignment_1.featureXML",maps[0]);
  xml_file.store("output/Tutorial_MapAlignment_2.featureXML",maps[1]);
  
  return 0;
} //end of main
