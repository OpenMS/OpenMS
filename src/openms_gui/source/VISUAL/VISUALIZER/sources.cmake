### the directory name
set(directory source/VISUAL/VISUALIZER)

### list all filenames of the directory here
set(sources_list
AcquisitionInfoVisualizer.cpp
AcquisitionVisualizer.cpp
BaseVisualizer.cpp
BaseVisualizerGUI.cpp
ContactPersonVisualizer.cpp
DataProcessingVisualizer.cpp
DocumentIdentifierVisualizer.cpp
ExperimentalSettingsVisualizer.cpp
GradientVisualizer.cpp
HPLCVisualizer.cpp
InstrumentSettingsVisualizer.cpp
InstrumentVisualizer.cpp
IonDetectorVisualizer.cpp
IonSourceVisualizer.cpp
MassAnalyzerVisualizer.cpp
MetaInfoDescriptionVisualizer.cpp
MetaInfoVisualizer.cpp
PeptideHitVisualizer.cpp
PeptideIdentificationVisualizer.cpp
PrecursorVisualizer.cpp
ProteinHitVisualizer.cpp
ProteinIdentificationVisualizer.cpp
SampleVisualizer.cpp
ScanWindowVisualizer.cpp
SoftwareVisualizer.cpp
SourceFileVisualizer.cpp
SpectrumSettingsVisualizer.cpp
ProductVisualizer.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL\\VISUALIZER" FILES ${sources})

