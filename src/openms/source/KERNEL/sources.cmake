### the directory name
set(directory source/KERNEL)

### list all filenames of the directory here
set(sources_list
AreaIterator.cpp
BaseFeature.cpp
ConsensusFeature.cpp
ConsensusMap.cpp
ConversionHelper.cpp
DPeak.cpp
Feature.cpp
FeatureHandle.cpp
FeatureMap.cpp
MassTrace.cpp
MRMFeature.cpp
MRMTransitionGroup.cpp
MSExperiment.cpp
MSSpectrum.cpp
OnDiscMSExperiment.cpp
Peak1D.cpp
Peak2D.cpp
PeakIndex.cpp
RangeManager.cpp
RichPeak2D.cpp
StandardTypes.cpp
ChromatogramPeak.cpp
MSChromatogram.cpp
ChromatogramTools.cpp
SpectrumHelper.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\KERNEL" FILES ${sources})

