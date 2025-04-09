### the directory name
set(directory include/OpenMS/ANALYSIS/MRM)

### list all header files of the directory here
set(sources_list_h
ReactionMonitoringTransition.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\MRM" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

