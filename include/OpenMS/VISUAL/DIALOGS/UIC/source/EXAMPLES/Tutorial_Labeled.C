#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

Int main()
{
  vector<FeatureMap<> > maps;
  maps.resize(1);

  FeatureXMLFile feature_file;
  feature_file.load("data/Tutorial_Labeled.featureXML",maps[0]);
  ConsensusMap out;
  out.getFileDescriptions()[0].filename = "data/Tutorial_Labeled.featureXML";
  out.getFileDescriptions()[0].size = maps[0].size();
  out.getFileDescriptions()[0].label = "light";
  out.getFileDescriptions()[1].filename = "data/Tutorial_Labeled.featureXML";
  out.getFileDescriptions()[1].size = maps[0].size();
  out.getFileDescriptions()[1].label = "heavy";

  FeatureGroupingAlgorithmLabeled algorithm;
  // ... set parameters
  algorithm.group(maps, out);
  ConsensusXMLFile consensus_file;
  consensus_file.store("output/Tutorial_Labeled.consensusXML",out);

  return 0;
} //end of main
