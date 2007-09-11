#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

using namespace OpenMS;
using namespace std;

Int main()
{
	FeatureFinder ff;
	
	Param parameters;
	// ... set parameters (e.g. from INI file)
	MSExperiment<> input;
	// ... set input data (e.g. from mzData file)
	FeatureMap<> output;

	ff.run("simple", input, output, parameters);
	
  return 0;
} //end of main
