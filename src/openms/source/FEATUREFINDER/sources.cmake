### the directory name
set(directory source/FEATUREFINDER)

### list all filenames of the directory here
set(sources_list
BaseModel.cpp
BiGaussFitter1D.cpp
BiGaussModel.cpp
EGHTraceFitter.cpp
ElutionModelFitter.cpp
ElutionPeakDetection.cpp
EmgFitter1D.cpp
EmgModel.cpp
EmgScoring.cpp
ExtendedIsotopeFitter1D.cpp
ExtendedIsotopeModel.cpp
FeatureFinderAlgorithmPicked.cpp
FeatureFinderAlgorithmPickedHelperStructs.cpp
FeatureFinderIdentificationAlgorithm.cpp
FeatureFinderAlgorithmMetaboIdent.cpp
FeatureFinderMultiplexAlgorithm.cpp
FeatureFindingMetabo.cpp
Fitter1D.cpp
GaussFitter1D.cpp
GaussModel.cpp
GaussTraceFitter.cpp
InterpolationModel.cpp
IsotopeFitter1D.cpp
IsotopeModel.cpp
LevMarqFitter1D.cpp
MassTraceDetection.cpp
MaxLikeliFitter1D.cpp
ModelDescription.cpp
MultiplexClustering.cpp
MultiplexDeltaMasses.cpp
MultiplexDeltaMassesGenerator.cpp
MultiplexFilteredMSExperiment.cpp
MultiplexFilteredPeak.cpp
MultiplexFiltering.cpp
MultiplexFilteringCentroided.cpp
MultiplexFilteringProfile.cpp
MultiplexIsotopicPeakPattern.cpp
MultiplexSatelliteCentroided.cpp
MultiplexSatelliteProfile.cpp
PeakWidthEstimator.cpp
SeedListGenerator.cpp
TraceFitter.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FEATUREFINDER" FILES ${sources})
