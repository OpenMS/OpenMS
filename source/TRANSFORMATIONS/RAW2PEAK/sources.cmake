### the directory name
set(directory source/TRANSFORMATIONS/RAW2PEAK)

### list all filenames of the directory here
set(sources_list
ContinuousWaveletTransform.C
ContinuousWaveletTransformNumIntegration.C
OptimizePeakDeconvolution.C
OptimizePick.C
PeakPickerCWT.C
PeakPickerHiRes.C
PeakShape.C
TwoDOptimization.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\TRANSFORMATIONS\\RAW2PEAK" FILES ${sources})

