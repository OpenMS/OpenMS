#include <OpenMS/KERNEL/FeatureMap.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  FeatureMap<> map;
  
  Feature feature;
  feature.setRT(15.0);
  feature.setMZ(571.3);
  map.push_back(feature); //append feature 1
  feature.setRT(23.3);
  feature.setMZ(1311.3);
  map.push_back(feature); //append feature 2
  
  
  for (FeatureMap<>::Iterator it=map.begin(); it!=map.end(); ++it)
  {
    cout << it->getRT() << " - " << it->getMZ() << endl;
  }

  return 0;
} //end of main
