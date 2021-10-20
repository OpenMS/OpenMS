### the directory name
set(directory source/FORMAT/OPTIONS)

### list all filenames of the directory here
set(sources_list
FeatureFileOptions.cpp
PeakFileOptions.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_FORMAT_sources ${OpenMS_FORMAT_sources} ${sources})

### source group definition
source_group("Source Files\\FORMAT\\OPTIONS" FILES ${sources})

