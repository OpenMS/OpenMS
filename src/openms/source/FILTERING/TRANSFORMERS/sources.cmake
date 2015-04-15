### the directory name
set(directory source/FILTERING/TRANSFORMERS)

### list all filenames of the directory here
set(sources_list
BernNorm.cpp
ComplementFilter.cpp
ComplementMarker.cpp
FilterFunctor.cpp
GoodDiffFilter.cpp
IntensityBalanceFilter.cpp
IsotopeDiffFilter.cpp
IsotopeMarker.cpp
LinearResampler.cpp
MarkerMower.cpp
NLargest.cpp
NeutralLossDiffFilter.cpp
NeutralLossMarker.cpp
Normalizer.cpp
ParentPeakMower.cpp
PeakMarker.cpp
#~ PreprocessingFunctor.cpp
Scaler.cpp
SpectraMerger.cpp
SqrtMower.cpp
TICFilter.cpp
ThresholdMower.cpp
WindowMower.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FILTERING\\TRANSFORMERS" FILES ${sources})

