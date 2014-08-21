### the directory name
set(directory source/TRANSFORMATIONS/FEATUREFINDER)

### list all filenames of the directory here
set(sources_list
BaseModel.cpp
BiGaussFitter1D.cpp
BiGaussModel.cpp
EGHTraceFitter.cpp
EmgFitter1D.cpp
EmgModel.cpp
EmgScoring.cpp
ExtendedIsotopeFitter1D.cpp
ExtendedIsotopeModel.cpp
FeaFiModule.cpp
FeatureFinder.cpp
FeatureFinderAlgorithm.cpp
FeatureFinderAlgorithmIsotopeWavelet.cpp
FeatureFinderAlgorithmMRM.cpp
FeatureFinderAlgorithmPicked.cpp
FeatureFinderAlgorithmPickedHelperStructs.cpp
FeatureFinderAlgorithmSimple.cpp
FeatureFinderAlgorithmSimplest.cpp
Fitter1D.cpp
GaussFitter1D.cpp
GaussModel.cpp
GaussTraceFitter.cpp
InterpolationModel.cpp
IsotopeFitter1D.cpp
IsotopeModel.cpp
IsotopeWavelet.cpp
IsotopeWaveletTransform.cpp
LevMarqFitter1D.cpp
MaxLikeliFitter1D.cpp
ModelDescription.cpp
ModelFitter.cpp
MultiplexClustering.cpp
MultiplexFiltering.cpp
MultiplexFilterResult.cpp
MultiplexFilterResultPeak.cpp
MultiplexFilterResultRaw.cpp
MultiplexPeakPattern.cpp
PeakWidthEstimator.cpp
ProductModel.cpp
SeedListGenerator.cpp
SimpleExtender.cpp
SimpleSeeder.cpp

)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\TRANSFORMATIONS\\FEATUREFINDER" FILES ${sources})
