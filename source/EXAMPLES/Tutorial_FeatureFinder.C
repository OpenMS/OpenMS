#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>

using namespace OpenMS;
using namespace std;

Int main()
{
  FeatureFinder ff;
  // ... set parameters (e.g. from INI file)
  Param parameters;
  // ... set input data (e.g. from mzML file)
  MSExperiment<> input;
  // ... set output data structure
  FeatureMap<> output;
  // ... set user-specified seeds, if needed
	FeatureMap<> seeds;

  ff.run("simple", input, output, parameters, seeds);
  
  return 0;
} //end of main
