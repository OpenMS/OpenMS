### the directory name
set(directory source/VISUAL/INTERFACES)

### list all filenames of the directory here
set(sources_list
  IPeptideIds.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL\\DIALOGS" FILES ${sources})

