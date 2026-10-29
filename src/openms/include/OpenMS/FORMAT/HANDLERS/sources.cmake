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
TraMLHandler.h
UnimodXMLHandler.h
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

### Private parser implementation and metadata headers, never installed or wrapped.
set(private_headers_list_h
ThermoRawFileMetadata.h
MzIdentMLDOMHandler.h
)
if (WITH_THERMO_RAW)
  list(APPEND private_headers_list_h ThermoRawScanReader.h)
endif()
set(private_sources_h)
foreach(i ${private_headers_list_h})
	list(APPEND private_sources_h ${directory}/${i})
endforeach(i)
source_group("Header Files\\OpenMS\\FORMAT\\HANDLERS" FILES ${private_sources_h})
set(OpenMS_private_headers ${OpenMS_private_headers} ${private_sources_h})

