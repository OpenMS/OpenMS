### the directory name
set(directory doc/doxygen/parameters)

### list all filenames of the directory here
set(sources_list
DefaultParamHandlerDocumenter
TOPPDocumenter
)

### pass source file list to the upper instance
set(OpenMS_doc_executables ${sources_list})

### source group definition
source_group("Source Files\\doc\\doxygen\\parameters" FILES ${sources})

