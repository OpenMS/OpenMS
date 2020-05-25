### the directory name
set(directory include/OpenMS/FORMAT/HANDLERS)

### list all header files of the directory here
set(sources_list_h
AcqusHandler.h
CachedMzMLHandler.h
FidHandler.h
IndexedMzMLDecoder.h
IndexedMzMLHandler.h
MascotXMLHandler.h
MzDataHandler.h
MzIdentMLDOMHandler.h
MzIdentMLHandler.h
MzMLHandler.h
MzMLHandlerHelper.h
MzMLSpectrumDecoder.h
MzMLSqliteHandler.h
MzMLSqliteSwathHandler.h
MzQuantMLHandler.h
MzXMLHandler.h
PTMXMLHandler.h
ParamXMLHandler.h
ToolDescriptionHandler.h
TraMLHandler.h
UnimodXMLHandler.h
XMLHandler.h
XQuestResultXMLHandler.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FORMAT\\HANDLERS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

