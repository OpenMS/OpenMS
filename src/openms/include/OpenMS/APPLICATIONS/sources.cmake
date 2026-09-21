### the directory name
set(directory include/OpenMS/APPLICATIONS)

### list all header files of the directory here
### (the TOPP tool framework -- TOPPBase, ToolHandler, ... -- lives in src/openms_cli)
set(sources_list_h
ConsoleUtils.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\APPLICATIONS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
