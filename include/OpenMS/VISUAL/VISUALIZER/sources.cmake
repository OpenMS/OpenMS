### the directory name
set(directory include/OpenMS/VISUAL/VISUALIZER)

### list all MOC filenames of the directory here
set(sources_list
AcquisitionInfoVisualizer.h
AcquisitionVisualizer.h
BaseVisualizerGUI.h
ContactPersonVisualizer.h
DataProcessingVisualizer.h
DigestionVisualizer.h
DocumentIdentifierVisualizer.h
ExperimentalSettingsVisualizer.h
GradientVisualizer.h
HPLCVisualizer.h
InstrumentSettingsVisualizer.h
InstrumentVisualizer.h
IonDetectorVisualizer.h
IonSourceVisualizer.h
MassAnalyzerVisualizer.h
MetaInfoDescriptionVisualizer.h
MetaInfoVisualizer.h
ModificationVisualizer.h
PeptideHitVisualizer.h
PeptideIdentificationVisualizer.h
PrecursorVisualizer.h
ProteinHitVisualizer.h
ProteinIdentificationVisualizer.h
SampleVisualizer.h
ScanWindowVisualizer.h
SoftwareVisualizer.h
SourceFileVisualizer.h
SpectrumSettingsVisualizer.h
TaggingVisualizer.h
ProductVisualizer.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply MOC compiler
QT4_WRAP_CPP(mocced_sources ${sources})

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${mocced_sources})

source_group("Source Files\\OpenMS\\VISUAL\\VISUALIZER" FILES ${mocced_sources})

### list all header files of the directory here
set(sources_list_h
AcquisitionInfoVisualizer.h
AcquisitionVisualizer.h
BaseVisualizer.h
BaseVisualizerGUI.h
ContactPersonVisualizer.h
DataProcessingVisualizer.h
DigestionVisualizer.h
DocumentIdentifierVisualizer.h
ExperimentalSettingsVisualizer.h
GradientVisualizer.h
HPLCVisualizer.h
InstrumentSettingsVisualizer.h
InstrumentVisualizer.h
IonDetectorVisualizer.h
IonSourceVisualizer.h
MassAnalyzerVisualizer.h
MetaInfoDescriptionVisualizer.h
MetaInfoVisualizer.h
ModificationVisualizer.h
PeptideHitVisualizer.h
PeptideIdentificationVisualizer.h
PrecursorVisualizer.h
ProteinHitVisualizer.h
ProteinIdentificationVisualizer.h
SampleVisualizer.h
ScanWindowVisualizer.h
SoftwareVisualizer.h
SourceFileVisualizer.h
SpectrumSettingsVisualizer.h
TaggingVisualizer.h
ProductVisualizer.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\VISUAL\\VISUALIZER" FILES ${sources_h})

set(OpenMSVisual_sources_h ${OpenMSVisual_sources_h} ${sources_h})

