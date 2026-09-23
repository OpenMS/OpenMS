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
CachedMzML.h
ChromeleonFile.h
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
FileInfo.h
FLASHDeconvFeatureFile.h
FLASHDeconvSpectrumFile.h
GNPSMetaValueFile.h
GNPSMGFFile.h
GNPSQuantificationFile.h
GzipIfstream.h
ZipIfstream.h
IBSpectraFile.h
IdXMLFile.h
ImzMLFile.h
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
OMSSACSVFile.h
OMSSAXMLFile.h
OSWFile.h
ParamCTDFile.h
ParamCWLFile.h
ParamJSONFile.h
ParamXMLFile.h
ParquetTableComparator.h
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
TransformationXMLFile.h
UnimodXMLFile.h
XMLFile.h
XTandemInfile.h
XTandemXMLFile.h
FileTypes.h
FileNameUtils.h
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
list(APPEND sources_list_h ArrowSchemaRegistry.h)
list(APPEND sources_list_h ArrowIOHelpers.h)
list(APPEND sources_list_h ParquetFilter.h)
list(APPEND sources_list_h XICParquetFile.h)
list(APPEND sources_list_h XIMParquetFile.h)
list(APPEND sources_list_h XIPMParquetFile.h)
list(APPEND sources_list_h QPXFile.h)
list(APPEND sources_list_h QPXIdentity.h)
list(APPEND sources_list_h ProteinGroupArrowExport.h)
list(APPEND sources_list_h QPXCollectionExport.h)
list(APPEND sources_list_h QPXValueValidation.h)
list(APPEND sources_list_h ProteinIdentificationArrowIO.h)
list(APPEND sources_list_h FeatureMapArrowIO.h)
list(APPEND sources_list_h ConsensusMapArrowIO.h)
list(APPEND sources_list_h PSMArrowIO.h)
list(APPEND sources_list_h ModificationDefinitionIO.h)

if (WITH_OPENTIMS)
  list(APPEND sources_list_h BrukerTimsFile.h)
  list(APPEND sources_list_h BrukerTimsImagingFile.h)
endif()

if (WITH_THERMO_RAW)
  list(APPEND sources_list_h ThermoRawFile.h)
endif()

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FORMAT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

### Private (non-installed) headers: the Xerces InputSource / BinInputStream
### adapters are internal plumbing used only by XMLFile.cpp / CompressedInputSource.cpp.
### Keeping them off OpenMS_sources_h is what lets Xerces be a PRIVATE link dependency.
###
### SqliteConnector_impl.h exposes the raw SQLite C API (sqlite3 / sqlite3_stmt)
### and OMSFileStore.h / OMSFileLoad.h expose the SQLiteCpp C++ API (SQLite::*).
### Keeping all three off OpenMS_sources_h is what lets SQLite (SQLiteCpp) be a
### fully private dependency: no SQLite type appears in any installed header.
###
### ParquetFile.h is the same case for Arrow: every one of its helpers takes or
### returns an arrow::Status / arrow::Table / arrow::Array, so the header includes
### <arrow/api.h> and cannot be compiled without Arrow's development files. It is
### an internal helper shared by the Parquet-backed I/O classes -- the installed
### readers/writers (XICParquetFile, QPXFile, ...) expose OpenMS types only -- so
### keeping it off OpenMS_sources_h is what lets Arrow/Parquet stay PRIVATE.
###
### ZipRandomAccessFile.h is the same: Open() returns an
### arrow::Result<std::shared_ptr<arrow::io::RandomAccessFile>>, so the header includes
### <arrow/io/api.h>. The Arrow-based OpenSWATH and Parquet helpers that use either one
### live in libOpenMS too (#10247), so no tool directory needs them installed.
set(private_headers_list_h
Bzip2InputStream.h
CompressedInputSource.h
GzipInputStream.h
ZipInputStream.h
SqliteConnector_impl.h
OMSFileLoad.h
OMSFileStore.h
ParquetFile.h
ZipRandomAccessFile.h
)

### RationalScan2ImConverter derives from OpenTIMS' Scan2InvIonMobilityConverter, so its
### header includes <opentims++/...> and needs OpenTIMS' development files. Only
### BrukerTimsFile.cpp and the class test use it; BrukerTimsFile.h itself hands out
### OpenMS types, so OpenTIMS stays PRIVATE.
if (WITH_OPENTIMS)
  list(APPEND private_headers_list_h RationalScan2ImConverter.h)
endif()

set(private_sources_h)
foreach(i ${private_headers_list_h})
	list(APPEND private_sources_h ${directory}/${i})
endforeach(i)
source_group("Header Files\\OpenMS\\FORMAT" FILES ${private_sources_h})
set(OpenMS_private_headers ${OpenMS_private_headers} ${private_sources_h})
