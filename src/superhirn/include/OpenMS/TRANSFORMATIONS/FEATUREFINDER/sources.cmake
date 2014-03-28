### the directory name
set(directory include/OpenMS/TRANSFORMATIONS/FEATUREFINDER)

### list all header files of the directory here
set(sources_list_h
FeatureFinderAlgorithmSH.h
FeatureFinderAlgorithmSHCtrl.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\TRANSFORMATIONS\\FEATUREFINDER" FILES ${sources_h})

set(SuperHirn_sources_h ${SuperHirn_sources_h} ${sources_h} CACHE INTERNAL "This variable should hold all SuperHirn headers at the end of the config step")
