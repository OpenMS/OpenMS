### the directory name
set(directory source/SIMULATION/LABELING)

### list all header files of the directory here
set(sources_list
BaseLabeler.cpp
ITRAQLabeler.cpp
LabelFreeLabeler.cpp
O18Labeler.cpp
SILACLabeler.cpp
ICPLLabeler.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_AUX_sources ${OpenMS_AUX_sources} ${sources})

### source group definition
source_group("Source Files\\SIMULATION\\LABELING" FILES ${sources})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)
