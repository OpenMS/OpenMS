### the directory name
set(directory source/APPLICATIONS)

### list all filenames of the directory here
set(sources_list
ConsoleUtils.cpp
INIUpdater.cpp
ParameterInformation.cpp
SearchEngineBase.cpp
ToolHandler.cpp
TOPPBase.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

if (MSVC)
	## apparently TOPPBase.cpp has more than 2^16 sections in its obj file (C1128). we need to increase to 2^32
	set_source_files_properties(${directory}/TOPPBase.cpp PROPERTIES COMPILE_FLAGS "/bigobj")
endif()

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\APPLICATIONS" FILES ${sources})

