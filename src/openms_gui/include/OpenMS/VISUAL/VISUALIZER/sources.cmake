### the directory name
set(directory include/OpenMS/VISUAL/VISUALIZER)

### list all header files of the directory here
set(sources_list_h
AcquisitionInfoVisualizer.h
AcquisitionVisualizer.h
BaseVisualizer.h               ### not a QObject -- should not be autoMOC'd
BaseVisualizerGUI.h
ContactPersonVisualizer.h
DataProcessingVisualizer.h
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
ProductVisualizer.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
  list(APPEND sources_h ${directory}/${i})
endforeach(i)

### treat as source files, for autoMOC'ing instead of manually calling QT5_WRAP_CPP()
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${sources_h})
### pass header file list to the upper instance
set(OpenMSVisual_sources_h ${OpenMSVisual_sources_h} ${sources_h})

### header group definition for IDE's
source_group("Header Files\\OpenMS\\VISUAL\\VISUALIZER" FILES ${sources_h})

