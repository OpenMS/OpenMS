### the directory name
set(directory source/ANALYSIS/TOPDOWN)

### list all filenames of the directory here
set(sources_list
MassFeatureTrace.cpp
DeconvolutedSpectrum.cpp
PeakGroup.cpp
QScore.cpp
PeakGroupScoring.cpp
SpectrumDeconvolution.cpp
FLASHDeconvAlgorithm.cpp
FLASHIda.cpp
FLASHDeconvHelperStructs.cpp
FLASHIdaBridgeFunctions.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\TOPDOWN" FILES ${sources})

