### the directory name
set(directory include/OpenMS/ML)

### list all filenames of the directory here
set(sources_list_h
	AlphapeptdeepWrapper.h AlphaDatahandling.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### pass source file list to the upper cmake instance
set(OpenMSML_sources_h ${OpenMSML_sources_h} ${sources_h})

### source group definition
source_group("Header Files\\OpenMSML" FILES ${sources_h})
