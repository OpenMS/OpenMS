### the directory name
set(directory source/FORMAT)

### list all filenames of the directory here
set(sources_list
AutoExecuteFile.C
Base64.C
CVMappingFile.C
CVMappings.C
ConsensusXMLFile.C
ControlledVocabulary.C
DTA2DFile.C
DTAFile.C
FASTAFile.C
FastaIterator.C
FastaIteratorIntern.C
FeatureXMLFile.C
FileHandler.C
IdXMLFile.C
InspectInfile.C
InspectOutfile.C
LibSVMEncoder.C
MS2File.C
MSPFile.C
MascotInfile.C
MascotInfile2.C
MascotOutfile.C
MascotRemoteQuery.C
MascotXMLFile.C
MzDataFile.C
MzMLFile.C
MzXMLFile.C
OMSSACSVFile.C
OMSSAXMLFile.C
PTMXMLFile.C
PeakFileOptions.C
PeakTypeEstimator.C
PepNovoInfile.C
PepNovoOutfile.C
PepXMLFile.C
ProtXMLFile.C
SequestInfile.C
SequestOutfile.C
TextFile.C
TransformationXMLFile.C
UnimodXMLFile.C
XMassFile.C
XMLFile.C
XTandemInfile.C
XTandemXMLFile.C
MzIdentMLFile.C
TraMLFile.C
)

if (USE_ANDIMS) 	 
	list(APPEND sources_list ANDIFile.C) 	 
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

