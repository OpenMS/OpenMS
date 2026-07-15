### the directory name
set(directory include/OpenMS/FORMAT/HANDLERS)

### list all header files of the directory here
set(sources_list_h
AcqusHandler.h
CachedMzMLHandler.h
ConsensusXMLHandler.h
FeatureXMLHandler.h
FidHandler.h
IndexedMzMLDecoder.h
ImzMLHandler.h
ImzMLHandlerHelper.h
ImzMLWriter.h
IndexedMzMLHandler.h
MascotXMLHandler.h
MzDataHandler.h
MzIdentMLHandler.h
MzMLHandler.h
MzMLHandlerHelper.h
MzMLSpectrumDecoder.h
MzMLSqliteHandler.h
MzMLSqliteSwathHandler.h
MzXMLHandler.h
PASEFHillCentroider.h
PTMXMLHandler.h
ParamXMLHandler.h
ToolDescriptionHandler.h
TraMLHandler.h
UnimodXMLHandler.h
UniProtXMLHandler.h
StringManager.h
XMLAttributes.h
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

### Private (non-installed) headers: internal to libOpenMS, never installed,
### exported, or wrapped. SAX2HandlerAdapter bridges Xerces SAX to the Xerces-free
### XMLHandler; MzIdentMLDOMHandler embeds a Xerces DOM parser (used only by
### MzIdentMLFile.cpp). Keeping them off OpenMS_sources_h is what lets Xerces be a
### PRIVATE link dependency.
set(private_headers_list_h
SAX2HandlerAdapter.h
MzIdentMLDOMHandler.h
)
set(private_sources_h)
foreach(i ${private_headers_list_h})
	list(APPEND private_sources_h ${directory}/${i})
endforeach(i)
source_group("Header Files\\OpenMS\\FORMAT\\HANDLERS" FILES ${private_sources_h})
set(OpenMS_private_headers ${OpenMS_private_headers} ${private_sources_h})

