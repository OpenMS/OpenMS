### the directory name
set(directory include/OpenMS/TRANSFORMATIONS/FEATUREFINDER)

### list all header files of the directory here
set(sources_list_h
BaseModel.h
BaseModel_impl.h
BiGaussFitter1D.h
BiGaussModel.h
EmgFitter1D.h
EmgModel.h
ExtendedIsotopeFitter1D.h
ExtendedIsotopeModel.h
FeaFiModule.h
FeatureFinder.h
FeatureFinderAlgorithm.h
FeatureFinderAlgorithmIsotopeWavelet.h
FeatureFinderAlgorithmPicked.h
FeatureFinderAlgorithmSimple.h
FeatureFinderAlgorithmSimplest.h
FeatureFinderAlgorithmWavelet.h
FeatureFinderAlgorithmMRM.h
FeatureFinderAlgorithm_impl.h
FeatureFinderDefs.h
FeatureFinder_impl.h
Fitter1D.h
GaussFitter1D.h
GaussModel.h
InterpolationModel.h
IsotopeFitter1D.h
IsotopeModel.h
IsotopeWavelet.h
IsotopeWaveletTransform.h
LevMarqFitter1D.h
LmaGaussFitter1D.h
LmaGaussModel.h
LmaIsotopeFitter1D.h
LmaIsotopeModel.h
MaxLikeliFitter1D.h
ModelDescription.h
ModelFitter.h
ProductModel.h
SimpleExtender.h
SimpleSeeder.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\TRANSFORMATIONS\\FEATUREFINDER" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

