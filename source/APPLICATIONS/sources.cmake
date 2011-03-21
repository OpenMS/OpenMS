### the directory name
set(directory source/APPLICATIONS)

### list all filenames of the directory here
set(sources_list
ToolHandler.C
TOPPBase.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

if (MSVC)
	## apparently TOPPBase.C has more than 2^16 sections in its obj file (C1128). we need to increase to 2^32
	SET_SOURCE_FILES_PROPERTIES(${directory}/TOPPBase.C	PROPERTIES COMPILE_FLAGS "/bigobj")
endif()

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\APPLICATIONS" FILES ${sources})

