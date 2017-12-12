### the directory name
set(directory source/FORMAT/HANDLERS)

### list all filenames of the directory here
set(sources_list
	AcqusHandler.cpp
	FidHandler.cpp
  IndexedMzMLDecoder.cpp
	MascotXMLHandler.cpp
	MzDataHandler.cpp
	MzIdentMLHandler.cpp
	MzIdentMLDOMHandler.cpp
	MzQuantMLHandler.cpp
	MzMLHandler.cpp
  MzMLHandlerHelper.cpp
  MzMLSpectrumDecoder.cpp
  MzMLSqliteHandler.cpp
  MzMLSqliteSwathHandler.cpp
	MzXMLHandler.cpp
	PTMXMLHandler.cpp
	ParamXMLHandler.cpp
	ToolDescriptionHandler.cpp
	TraMLHandler.cpp
	UnimodXMLHandler.cpp
	XMLHandler.cpp
	XQuestResultXMLHandler.cpp
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

