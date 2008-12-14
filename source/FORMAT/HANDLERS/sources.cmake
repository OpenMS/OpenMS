### the directory name
set(directory source/FORMAT/HANDLERS)

### list all filenames of the directory here
set(sources_list
ANDIHandler.C
MascotXMLHandler.C
MzDataHandler.C
MzMLHandler.C
MzXMLHandler.C
PTMXMLHandler.C
ParamXMLHandler.C
UnimodXMLHandler.C
XMLHandler.C
XTandemInfileXMLHandler.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FORMAT\\HANDLERS" FILES ${sources})

