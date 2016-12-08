### the directory name
set(directory source/TRANSFORMATIONS/FEATUREFINDER)

### list all filenames of the directory here
set(sources_list
FeatureFinderAlgorithmSH.cpp
FeatureFinderAlgorithmSHCtrl.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(SuperHirn_sources ${SuperHirn_sources} ${sources} CACHE INTERNAL "This variable should hold all SuperHirn sources at the end of the config step" )

### source group definition
source_group("Source Files\\TRANSFORMATIONS\\FEATUREFINDER" FILES ${sources})
  
