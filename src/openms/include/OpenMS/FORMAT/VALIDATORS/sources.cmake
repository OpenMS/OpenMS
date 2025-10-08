### the directory name
set(directory include/OpenMS/FORMAT/VALIDATORS)

### list all header files of the directory here
set(sources_list_h
XMLValidator.h
MzMLValidator.h
MzDataValidator.h
SemanticValidator.h
MzIdentMLValidator.h
TraMLValidator.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FORMAT\\VALIDATORS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

