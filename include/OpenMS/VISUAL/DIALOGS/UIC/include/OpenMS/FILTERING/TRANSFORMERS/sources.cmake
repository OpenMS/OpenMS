### the directory name
set(directory include/OpenMS/FILTERING/TRANSFORMERS)

### list all header files of the directory here
set(sources_list_h
BernNorm.h
ComplementFilter.h
ComplementMarker.h
FilterFunctor.h
GoodDiffFilter.h
IntensityBalanceFilter.h
IsotopeDiffFilter.h
IsotopeMarker.h
LinearResampler.h
MarkerMower.h
NLargest.h
NeutralLossDiffFilter.h
NeutralLossMarker.h
Normalizer.h
ParentPeakMower.h
PeakMarker.h
PreprocessingFunctor.h
Scaler.h
SqrtMower.h
TICFilter.h
ThresholdMower.h
WindowMower.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FILTERING\\TRANSFORMERS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

