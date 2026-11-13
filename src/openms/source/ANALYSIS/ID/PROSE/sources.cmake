### the directory name
set(directory source/ANALYSIS/ID/PROSE)

### list all filenames of the directory here
set(sources_list
ProSEAlgorithm.cpp
ProSEHits.cpp
ProSELoading.cpp
ProSEPEFF.cpp
ProSEPooling.cpp
ProSEProteins.cpp
ProSEReport.cpp
ProSERescoring.cpp
ProSEScoring.cpp
ProSESpectra.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### private (non-installed) headers of this directory
list(APPEND OpenMS_private_headers ${directory}/ProSEImpl.h)

### source group definition
source_group("Source Files\\ANALYSIS\\ID\\PROSE" FILES ${sources})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)
