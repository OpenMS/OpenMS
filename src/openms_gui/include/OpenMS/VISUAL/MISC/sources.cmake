### the directory name
set(directory include/OpenMS/VISUAL/MISC)

### list all header files of the directory here
set(sources_list_h
ExternalProcessMBox.h
GUIHelpers.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\VISUAL\\MISC" FILES ${sources_h})

set(OpenMSVisual_sources_h ${OpenMSVisual_sources_h} ${sources_h})

