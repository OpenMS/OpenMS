### the directory name
set(directory include/OpenMS/OPENSWATHALGO/DATAACCESS)

### list all header files of the directory here
set(sources_list_h
DataFrameWriter.h
DataStructures.h
ISpectrumAccess.h
ITrans2Trans.h
ITransition.h
MockObjects.h
TransitionExperiment.h
Transitions.h
SpectrumHelpers.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\OPENSWATH\\OPENSWATHALGO\\DATAACCESS" FILES ${sources_h})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)

set(OpenSWATHALGO_DATAACCESS_sources_h ${sources_h})
set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

