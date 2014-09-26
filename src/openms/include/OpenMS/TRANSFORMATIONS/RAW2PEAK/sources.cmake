### the directory name
set(directory include/OpenMS/TRANSFORMATIONS/RAW2PEAK)

### list all MOC filenames of the directory here
set(sources_list
PeakInvestigator.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply MOC compiler
QT4_WRAP_CPP(mocced_sources ${sources} OPTIONS ${BOOST_MOC_ARGS})

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${mocced_sources})

source_group("Source Files\\OpenMS\\TRANSFORMATIONS\\RAW2PEAK" FILES ${mocced_sources})

### list all header files of the directory here
set(sources_list_h
ContinuousWaveletTransform.h
ContinuousWaveletTransformNumIntegration.h
OptimizePeakDeconvolution.h
OptimizePick.h
PeakInvestigator.h
PeakPickerCWT.h
PeakPickerHiRes.h
PeakPickerIterative.h
PeakPickerSH.h
PeakShape.h
TwoDOptimization.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\TRANSFORMATIONS\\RAW2PEAK" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

