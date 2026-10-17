### the directory name
set(directory include/OpenMS/VISUAL/VISITORS)

### list all header files of the directory here
set(sources_list_h
LayerStatistics.h
LayerStoreData.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### pass header file list to the upper instance
set(OpenMSVisual_sources_h ${OpenMSVisual_sources_h} ${sources_h})

### header group definition for IDE's
source_group("Header Files\\OpenMS\\VISUAL\\VISITORS" FILES ${sources_h})
