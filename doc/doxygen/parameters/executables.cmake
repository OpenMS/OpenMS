### the directory name
set(directory doc/doxygen/parameters)

### list all filenames of the directory here
set(sources_list
DefaultParamHandlerDocumenter
TOPPDocumenter
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_doc_executables ${sources})

### source group definition
source_group("Source Files\\doc\\doxygen\\parameters" FILES ${sources})

