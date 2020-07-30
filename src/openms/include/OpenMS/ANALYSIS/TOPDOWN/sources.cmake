### the directory name
set(directory include/OpenMS/ANALYSIS/TOPDOWN)

### list all header files of the directory here
set(sources_list_h
MassFeatureTrace.h
DeconvolutedSpectrum.h
PeakGroup.h
		QScore.h
PeakGroupScoring.h
SpectrumDeconvolution.h
FLASHDeconvAlgorithm.h
FLASHIda.h
FLASHDeconvHelperStructs.h
FLASHIdaBridgeFunctions.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\TOPDOWN" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

