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
FASTAFile.cpp
FastaIterator.cpp
FastaIteratorIntern.cpp
FeatureXMLFile.cpp
FileHandler.cpp
FileTypes.cpp
GzipIfstream.cpp
GzipInputStream.cpp
IBSpectraFile.cpp
IdXMLFile.cpp
IndexedMzMLFile.cpp
IndexedMzMLFileLoader.cpp
InspectInfile.cpp
InspectOutfile.cpp
KroenikFile.cpp
LibSVMEncoder.cpp
MRMFeatureQCFile.cpp
MS2File.cpp
MSNumpressCoder.cpp
MSPFile.cpp
MascotInfile.cpp
MascotGenericFile.cpp
MascotRemoteQuery.cpp
MascotXMLFile.cpp
MsInspectFile.cpp
MzDataFile.cpp
MzIdentMLFile.cpp
MzMLFile.cpp
MzQuantMLFile.cpp
MzTab.cpp
MzTabFile.cpp
MzXMLFile.cpp
OMSSACSVFile.cpp
OMSSAXMLFile.cpp
OSWFile.cpp
ParamXMLFile.cpp
PTMXMLFile.cpp
PeakTypeEstimator.cpp
PepNovoInfile.cpp
PepNovoOutfile.cpp
PepXMLFile.cpp
PepXMLFileMascot.cpp
PercolatorOutfile.cpp
ProtXMLFile.cpp
QcMLFile.cpp
SequestInfile.cpp
SequestOutfile.cpp
SpecArrayFile.cpp
SqMassFile.cpp
SwathFile.cpp
SVOutStream.cpp
TextFile.cpp
ToolDescriptionFile.cpp
TraMLFile.cpp
TransformationXMLFile.cpp
UnimodXMLFile.cpp
XMassFile.cpp
XMLFile.cpp
XQuestResultXMLFile.cpp
XTandemInfile.cpp
XTandemXMLFile.cpp
ZlibCompression.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FORMAT" FILES ${sources})

