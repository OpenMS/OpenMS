### the directory name
set(directory include/OpenMS/FORMAT)

### list all MOC filenames of the directory here
set(sources_list
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
FLASHDeconvFeatureFile.h
FLASHDeconvSpectrumFile.h
GNPSMetaValueFile.h
GNPSMGFFile.h
GNPSQuantificationFile.h
GzipIfstream.h
GzipInputStream.h
ZipIfstream.h
ZipInputStream.h
IBSpectraFile.h
IdXMLFile.h
IndentedStream.h
IndexedMzMLFileLoader.h
InspectInfile.h
InspectOutfile.h
KroenikFile.h
MRMFeaturePickerFile.h
MRMFeatureQCFile.h
MS2File.h
MSNumpressCoder.h
MSPFile.h
MSPGenericFile.h
MSstatsFile.h
MascotGenericFile.h
MascotRemoteQuery.h
MascotXMLFile.h
MsInspectFile.h
MzDataFile.h
MzMLFile.h
MzQCFile.h
MzTab.h
MzTabBase.h
MzTabM.h
MzTabFile.h
MzTabMFile.h
MzXMLFile.h
OMSFile.h
OMSFileLoad.h
OMSFileStore.h
OMSSACSVFile.h
OMSSAXMLFile.h
OSWFile.h
ParamCTDFile.h
ParamCWLFile.h
ParamJSONFile.h
ParamXMLFile.h
PEFFFile.h
PTMXMLFile.h
PeakTypeEstimator.h
PepNovoInfile.h
PepNovoOutfile.h
PepXMLFile.h
PepXMLFileMascot.h
PercolatorInfile.h
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
TriqlerFile.h
UnimodXMLFile.h
XMLFile.h
XTandemInfile.h
XTandemXMLFile.h
FileTypes.h
MzIdentMLFile.h
TraMLFile.h
XMassFile.h
XQuestResultXMLFile.h
MRMFile.h
TargetedDataFileLoader.h
ZlibCompression.h
)

if (WITH_HDF5)
  list(APPEND sources_list_h HDF5Connector.h)
endif()

list(APPEND sources_list_h ZipArchiveFile.h)
list(APPEND sources_list_h MSExperimentArrowExport.h)
list(APPEND sources_list_h ConsensusMapArrowExport.h)
list(APPEND sources_list_h ParquetFile.h)
list(APPEND sources_list_h ParquetFilter.h)
list(APPEND sources_list_h XICParquetFile.h)
list(APPEND sources_list_h XIMParquetFile.h)
list(APPEND sources_list_h QPXFile.h)
list(APPEND sources_list_h ProteinGroupArrowExport.h)
list(APPEND sources_list_h ProteinIdentificationArrowIO.h)
list(APPEND sources_list_h FeatureMapArrowIO.h)
list(APPEND sources_list_h ConsensusMapArrowIO.h)

if (WITH_OPENTIMS)
  list(APPEND sources_list_h BrukerTimsFile.h)
  list(APPEND sources_list_h RationalScan2ImConverter.h)
endif()

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FORMAT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
