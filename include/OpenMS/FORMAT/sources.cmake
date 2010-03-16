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

### Apply MOC compiler
QT4_WRAP_CPP(mocced_sources ${sources})

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${mocced_sources})

source_group("Source Files\\OpenMS\\FORMAT" FILES ${mocced_sources})

### list all header files of the directory here
set(sources_list_h
ANDIFile.h
Base64.h
Bzip2Ifstream.h
Bzip2InputStream.h
CompressedInputSource.h
CVMappingFile.h
ConsensusXMLFile.h
ControlledVocabulary.h
CsvFile.h
DTA2DFile.h
DTAFile.h
EDTAFile.h
FASTAFile.h
FastaIterator.h
FastaIteratorIntern.h
FeatureXMLFile.h
FileHandler.h
GzipIfstream.h
GzipInputStream.h
IdXMLFile.h
InspectInfile.h
InspectOutfile.h
KroenikFile.h
LibSVMEncoder.h
MS2File.h
MSPFile.h
MascotInfile.h
MascotGenericFile.h
MascotOutfile.h
MascotRemoteQuery.h
MascotXMLFile.h
MsInspectFile.h
MzDataFile.h
MzMLFile.h
MzXMLFile.h
OMSSACSVFile.h
OMSSAXMLFile.h
PTMXMLFile.h
PeakFileOptions.h
PeakTypeEstimator.h
PepNovoInfile.h
PepNovoOutfile.h
PepXMLFile.h
PepXMLFileMascot.h
ProtXMLFile.h
SequestInfile.h
SequestOutfile.h
SpecArrayFile.h
SVOutStream.h
TextFile.h
TransformationXMLFile.h
UnimodXMLFile.h
XMLFile.h
XTandemInfile.h
XTandemXMLFile.h
FileTypes.h
MzIdentMLFile.h
TraMLFile.h
XMassFile.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FORMAT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

