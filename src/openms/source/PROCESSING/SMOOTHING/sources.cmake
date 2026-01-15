### the directory name
set(directory source/PROCESSING/SMOOTHING)

### list all filenames of the directory here
set(sources_list
GaussFilter.cpp
GaussFilterAlgorithm.cpp
LowessSmoothing.cpp
FastLowessSmoothing.cpp
SavitzkyGolayFilter.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FILTERING\\SMOOTHING" FILES ${sources})

