### the directory name
set(directory source/FORMAT)

### list all filenames of the directory here
set(sources_list
AbsoluteQuantitationMethodFile.cpp
AbsoluteQuantitationStandardsFile.cpp
Base64.cpp
Bzip2Ifstream.cpp
Bzip2InputStream.cpp
CachedMzML.cpp
ChromeleonFile.cpp
CompressedInputSource.cpp
CVMappingFile.cpp
ConsensusXMLFile.cpp
ControlledVocabulary.cpp
CsvFile.cpp
DTA2DFile.cpp
DTAFile.cpp
EDTAFile.cpp
ExperimentalDesignFile.cpp
FASTAFile.cpp
FeatureXMLFile.cpp
FLASHDeconvFeatureFile.cpp
FLASHDeconvSpectrumFile.cpp
FileHandler.cpp
FileTypes.cpp
GNPSMetaValueFile.cpp
GNPSMGFFile.cpp
GNPSQuantificationFile.cpp
GzipIfstream.cpp
GzipInputStream.cpp
ZipIfstream.cpp
ZipInputStream.cpp
IBSpectraFile.cpp
IdXMLFile.cpp
IndentedStream.cpp
IndexedMzMLFileLoader.cpp
InspectInfile.cpp
InspectOutfile.cpp
KroenikFile.cpp
MRMFeaturePickerFile.cpp
MRMFeatureQCFile.cpp
MS2File.cpp
MSNumpressCoder.cpp
MSPFile.cpp
MSPGenericFile.cpp
MSstatsFile.cpp
MascotGenericFile.cpp
MascotRemoteQuery.cpp
MascotXMLFile.cpp
MsInspectFile.cpp
MzDataFile.cpp
MzIdentMLFile.cpp
MzMLFile.cpp
MzQCFile.cpp
MzTab.cpp
MzTabBase.cpp
MzTabM.cpp
MzTabFile.cpp
MzTabMFile.cpp
MzXMLFile.cpp
OMSFile.cpp
OMSFileLoad.cpp
OMSFileStore.cpp
OMSSACSVFile.cpp
OMSSAXMLFile.cpp
OSWFile.cpp
ParamCTDFile.cpp
ParamCWLFile.cpp
ParamJSONFile.cpp
ParamXMLFile.cpp
PEFFFile.cpp
PTMXMLFile.cpp
PeakTypeEstimator.cpp
PepNovoInfile.cpp
PepNovoOutfile.cpp
PepXMLFile.cpp
PepXMLFileMascot.cpp
PercolatorInfile.cpp
PercolatorOutfile.cpp
ProtXMLFile.cpp
QcMLFile.cpp
SequestInfile.cpp
SequestOutfile.cpp
SpecArrayFile.cpp
SqliteConnector.cpp
SqMassFile.cpp
SwathFile.cpp
SVOutStream.cpp
TextFile.cpp
ToolDescriptionFile.cpp
TraMLFile.cpp
TransformationXMLFile.cpp
TriqlerFile.cpp
UnimodXMLFile.cpp
XMassFile.cpp
XMLFile.cpp
XQuestResultXMLFile.cpp
XTandemInfile.cpp
XTandemXMLFile.cpp
ZlibCompression.cpp
MRMFile.cpp
TargetedDataFileLoader.cpp
)

if (WITH_HDF5)
  list(APPEND sources_list HDF5Connector.cpp)
endif()

list(APPEND sources_list ZipArchiveFile.cpp)
list(APPEND sources_list ZipRandomAccessFile.cpp)
list(APPEND sources_list MSExperimentArrowExport.cpp)
list(APPEND sources_list ConsensusMapArrowExport.cpp)
list(APPEND sources_list ParquetFile.cpp)
list(APPEND sources_list ParquetFilter.cpp)
list(APPEND sources_list QPXFile.cpp)
list(APPEND sources_list XICParquetFile.cpp)
list(APPEND sources_list XIMParquetFile.cpp)
list(APPEND sources_list ProteinGroupArrowExport.cpp)
list(APPEND sources_list ProteinIdentificationArrowIO.cpp)
list(APPEND sources_list FeatureMapArrowIO.cpp)
list(APPEND sources_list ConsensusMapArrowIO.cpp)
list(APPEND sources_list PSMArrowIO.cpp)
list(APPEND sources_list ArrowSchemaRegistry.cpp)
list(APPEND sources_list ArrowIOHelpers.cpp)

if (WITH_OPENTIMS)
  list(APPEND sources_list BrukerTimsFile.cpp)
  list(APPEND sources_list BrukerTimsImagingFile.cpp)
  list(APPEND sources_list RationalScan2ImConverter.cpp)
endif()

if (WITH_THERMO_RAW)
  list(APPEND sources_list ThermoRawFile.cpp)
endif()

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FORMAT" FILES ${sources})
