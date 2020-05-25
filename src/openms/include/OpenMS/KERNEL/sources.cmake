### the directory name
set(directory include/OpenMS/KERNEL)

### list all header files of the directory here
set(sources_list_h
AreaIterator.h
BaseFeature.h
ChromatogramPeak.h
ChromatogramTools.h
ComparatorUtils.h
ConsensusFeature.h
ConversionHelper.h
ConsensusMap.h
ConversionHelper.h
DPeak.h
Feature.h
FeatureHandle.h
FeatureMap.h
MassTrace.h
MRMFeature.h
MRMTransitionGroup.h
MSChromatogram.h
MSExperiment.h
MSSpectrum.h
OnDiscMSExperiment.h
Peak1D.h
Peak2D.h
PeakIndex.h
RangeManager.h
RangeUtils.h
RichPeak2D.h
StandardTypes.h
StandardDeclarations.h
SpectrumHelper.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\KERNEL" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

