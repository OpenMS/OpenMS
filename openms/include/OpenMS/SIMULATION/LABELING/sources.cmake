### the directory name
set(directory include/OpenMS/SIMULATION/LABELING)

### list all header files of the directory here
set(sources_list_h
BaseLabeler.h
BaseLabeler_impl.h
ITRAQLabeler.h
LabelFreeLabeler.h
O18Labeler.h
SILACLabeler.h
ICPLLabeler.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
  list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\SIMULATION\\LABELING" FILES ${sources_h})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

