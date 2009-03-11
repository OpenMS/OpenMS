#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

Int main()
{
  vector<FeatureMap<> > maps;
  maps.resize(2);

  FeatureXMLFile feature_file;
  feature_file.load("data/Tutorial_Unlabeled_1.featureXML",maps[0]);
  feature_file.load("data/Tutorial_Unlabeled_2.featureXML",maps[1]);

  ConsensusMap out;
  out.getFileDescriptions()[0].filename = "data/Tutorial_Unlabeled_1.featureXML";
  out.getFileDescriptions()[0].size = maps[0].size();
  out.getFileDescriptions()[1].filename = "data/Tutorial_Unlabeled_2.featureXML";
  out.getFileDescriptions()[1].size = maps[1].size();

  FeatureGroupingAlgorithmUnlabeled algorithm;
  // ... set parameters
  algorithm.group(maps, out);
  ConsensusXMLFile consensus_file;
  consensus_file.store("output/Tutorial_Unlabeled.consensusXML",out);

  return 0;
} //end of main
