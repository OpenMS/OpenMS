### the directory name
set(directory include/OpenMS/FORMAT)

### list all MOC filenames of the directory here
set(sources_list
MascotRemoteQuery.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### treat as source files, for autoMOC'ing instead of manually calling QT5_WRAP_CPP()
set(OpenMS_sources ${OpenMS_sources} ${sources})

### list all header files of the directory here
set(sources_list_h
AbsoluteQuantitationMethodFile.h
AbsoluteQuantitationStandardsFile.h
Base64.h
Bzip2Ifstream.h
Bzip2InputStream.h
CachedMzML.h
ChromeleonFile.h
CompressedInputSource.h
CVMappingFile.h
ConsensusXMLFile.h
ControlledVocabulary.h
CsvFile.h
DTA2DFile.h
DTAFile.h
EDTAFile.h
ExperimentalDesignFile.h
FASTAFile.h
FeatureXMLFile.h
FileHandler.h
GzipIfstream.h
GzipInputStream.h
HDF5Connector.h
IBSpectraFile.h
IdXMLFile.h
IndexedMzMLFileLoader.h
InspectInfile.h
InspectOutfile.h
KroenikFile.h
LibSVMEncoder.h
MRMFeaturePickerFile.h
MRMFeatureQCFile.h
MS2File.h
MSNumpressCoder.h
MSPFile.h
MSPGenericFile.h
MSstatsFile.h
MascotInfile.h
MascotGenericFile.h
MascotRemoteQuery.h
MascotXMLFile.h
MsInspectFile.h
MzDataFile.h
MzMLFile.h
MzTab.h
MzTabFile.h
MzXMLFile.h
OMSSACSVFile.h
OMSSAXMLFile.h
OSWFile.h
ParamXMLFile.h
PTMXMLFile.h
PeakTypeEstimator.h
PepNovoInfile.h
PepNovoOutfile.h
PepXMLFile.h
PepXMLFileMascot.h
PercolatorOutfile.h
ProtXMLFile.h
QcMLFile.h
SequestInfile.h
SequestOutfile.h
SpecArrayFile.h
SVOutStream.h
SwathFile.h
SqliteConnector.h
SqMassFile.h
TextFile.h
ToolDescriptionFile.h
TransformationXMLFile.h
UnimodXMLFile.h
XMLFile.h
XTandemInfile.h
XTandemXMLFile.h
FileTypes.h
MzIdentMLFile.h
MzQuantMLFile.h
TraMLFile.h
XMassFile.h
XQuestResultXMLFile.h
ZlibCompression.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FORMAT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
