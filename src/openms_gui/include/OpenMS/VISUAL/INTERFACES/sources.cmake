### the directory name
set(directory include/OpenMS/VISUAL/INTERFACES)

### list all header files of the directory here
set(sources_list_h
  IPeptideIds.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### pass header file list to the upper instance
set(OpenMSVisual_sources_h ${OpenMSVisual_sources_h} ${sources_h})

### header group definition for IDE's
source_group("Header Files\\OpenMS\\VISUAL\\INTERFACES" FILES ${sources_h})
