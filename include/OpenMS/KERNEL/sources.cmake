### the directory name
set(directory include/OpenMS/KERNEL)

### list all header files of the directory here
set(sources_list_h
AreaIterator.h
BaseFeature.h
ComparatorUtils.h
ConsensusFeature.h
ConsensusMap.h
DPeak.h
DRichPeak.h
Feature.h
FeatureHandle.h
FeatureMap.h
MSExperiment.h
MSSpectrum.h
Peak1D.h
Peak2D.h
PeakIndex.h
RangeManager.h
RangeUtils.h
RichPeak1D.h
RichPeak2D.h
StandardTypes.h
ChromatogramPeak.h
MSChromatogram.h
ChromatogramTools.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\KERNEL" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

