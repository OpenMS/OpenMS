### the directory name
set(directory include/OpenMS/TRANSFORMATIONS/FEATUREFINDER)

### list all header files of the directory here
set(sources_list_h
BaseModel.h
BaseModel_impl.h
BiGaussFitter1D.h
BiGaussModel.h
EGHTraceFitter.h
ElutionModelFitter.h
EmgFitter1D.h
EmgModel.h
EmgScoring.h
ExtendedIsotopeFitter1D.h
ExtendedIsotopeModel.h
FeatureFinder.h
FeatureFinderAlgorithm.h
FeatureFinderAlgorithmIsotopeWavelet.h
FeatureFinderAlgorithmMRM.h
FeatureFinderAlgorithmPicked.h
FeatureFinderAlgorithmPickedHelperStructs.h
FeatureFinderIdentificationAlgorithm.h
FeatureFinderMultiplexAlgorithm.h
FeatureFinderDefs.h
Fitter1D.h
GaussFitter1D.h
GaussModel.h
GaussTraceFitter.h
InterpolationModel.h
IsotopeFitter1D.h
IsotopeModel.h
IsotopeWavelet.h
IsotopeWaveletConstants.h
IsotopeWaveletTransform.h
LevMarqFitter1D.h
MaxLikeliFitter1D.h
ModelDescription.h
MultiplexClustering.h
MultiplexDeltaMasses.h
MultiplexDeltaMassesGenerator.h
MultiplexFilteredMSExperiment.h
MultiplexFilteredPeak.h
MultiplexFiltering.h
MultiplexFilteringCentroided.h
MultiplexFilteringProfile.h
MultiplexIsotopicPeakPattern.h
MultiplexSatelliteCentroided.h
MultiplexSatelliteProfile.h
PeakWidthEstimator.h
ProductModel.h
SeedListGenerator.h
TraceFitter.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
  list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\TRANSFORMATIONS\\FEATUREFINDER" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
