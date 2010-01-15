### the directory name
set(directory source/TRANSFORMATIONS/FEATUREFINDER)

### list all filenames of the directory here
set(sources_list
BaseModel.C
BiGaussFitter1D.C
BiGaussModel.C
EmgFitter1D.C
EmgModel.C
ExtendedIsotopeFitter1D.C
ExtendedIsotopeModel.C
FeaFiModule.C
FeatureFinder.C
FeatureFinderAlgorithm.C
FeatureFinderAlgorithmIsotopeWavelet.C
FeatureFinderAlgorithmPicked.C
FeatureFinderAlgorithmSimple.C
FeatureFinderAlgorithmSimplest.C
FeatureFinderAlgorithmMRM.C
Fitter1D.C
GaussFitter1D.C
GaussModel.C
InterpolationModel.C
IsotopeFitter1D.C
IsotopeModel.C
IsotopeWavelet.C
LevMarqFitter1D.C
LmaGaussFitter1D.C
LmaGaussModel.C
LmaIsotopeFitter1D.C
LmaIsotopeModel.C
MaxLikeliFitter1D.C
ModelDescription.C
ModelFitter.C
ProductModel.C
SeedListGenerator.C
SimpleExtender.C
SimpleSeeder.C
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

set(sources_list_cu
IsotopeWaveletCudaKernel.cu
)

set(sources_cu)
foreach(i ${sources_list_cu})
        list(APPEND sources_cu ${directory}/${i})
endforeach(i)

### 
set(Cuda_sources ${Cuda_sources} ${sources_cu})
