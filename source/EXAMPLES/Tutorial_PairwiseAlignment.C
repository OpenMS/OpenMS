#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringPairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapMatcherRegression.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
// #include <OpenMS/FORMAT/GridFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  FeatureMap<> exp_feature_1;
  FeatureMap<> exp_feature_2;
  
  FeatureXMLFile featurexml_file;
  featurexml_file.load("../TEST/TOPP/MapAlignmentFeatureMap1.xml",exp_feature_1);
  featurexml_file.load("../TEST/TOPP/MapAlignmentFeatureMap2.xml",exp_feature_2);
  
  Param param;
  param.setValue("superimposer:type","poseclustering_affine");
  param.setValue("superimposer:tuple_search:mz_bucket_size",0.3);
  param.setValue("pairfinder:type","DelaunayPairFinder");
  
  std::vector < ElementPair<Feature> > landmarks;
  PoseClusteringPairwiseMapMatcher< FeatureMap<> > pcpm;
  pcpm.setParameters(param);
  pcpm.setElementMap(0,exp_feature_1);    
  pcpm.setElementMap(1,exp_feature_2);
  pcpm.run();
  
	// TODO adapt to new map alignment class structure
#if 0
	GridFile grid_file;
  grid_file.store("FirstAffineTransformation.gridXML",pcpm.getGrid());
  
  MapMatcherRegression< Feature > lr;
  lr.setElementPairs(pcpm.getElementPairs());
  lr.setGrid(pcpm.getGrid());
  lr.estimateTransform();
  grid_file.store("SecondAffineTransformation.gridXML",lr.getGrid());    
#endif

  return 0;
} //end of main
