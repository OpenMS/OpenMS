#include <OpenMS/KERNEL/FeatureMap.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  FeatureMap<> map;

  Feature feature;
  feature.setIntensity(461.3f);
  feature.setRT(15.0);
  feature.setMZ(571.3);
  map.push_back(feature);
  feature.setIntensity(12213.5f);
  feature.setRT(23.3);
  feature.setMZ(1311.3);
  map.push_back(feature);

  //calculate the ranges
  map.updateRanges();

  cout << "Int: " << map.getMinInt() << " - " << map.getMaxInt() << endl;
  cout << "RT:  " << map.getMin()[0] << " - " << map.getMax()[0] << endl;
  cout << "m/z: " << map.getMin()[1] << " - " << map.getMax()[1] << endl;

  return 0;
} //end of main
