### the directory name
set(directory source/ML)

### list all filenames of the directory here
set(sources_list
	OpenMSMLWrapper.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMSML_sources ${OpenMSML_sources} ${sources})

### source group definition
source_group("Source Files\\OpenMSML" FILES ${sources})
