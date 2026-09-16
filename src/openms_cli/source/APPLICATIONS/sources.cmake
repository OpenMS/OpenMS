### the directory name
set(directory source/APPLICATIONS)

### list all filenames of the directory here
set(sources_list
INIUpdater.cpp
MapAlignerBase.cpp
OpenSwathBase.cpp
ParameterInformation.cpp
SearchEngineBase.cpp
ToolHandler.cpp
TOPPBase.cpp
TOPPExternalToolBase.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_CLI_sources ${OpenMS_CLI_sources} ${sources})

### source group definition
source_group("Source Files\\APPLICATIONS" FILES ${sources})
