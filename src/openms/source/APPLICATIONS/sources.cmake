### the directory name
set(directory source/APPLICATIONS)

### list all filenames of the directory here
### (the TOPP tool framework -- TOPPBase, ToolHandler, ... -- lives in src/openms_cli)
set(sources_list
ConsoleUtils.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\APPLICATIONS" FILES ${sources})
