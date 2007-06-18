#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

using namespace OpenMS;
using namespace std;

Int main()
{
	FeatureFinder ff;
	ff.addSeeder("SimpleSeeder");
	ff.addExtender("SimpleExtender");
	ff.addFitter("SimpleModelFitter");

	MSExperiment<> exp;
	//... Fill peak map with points
	ff.setData(exp.begin(),exp.end(),300);
	
	FeatureMap<> ouput = ff.run();

  return 0;
} //end of main
