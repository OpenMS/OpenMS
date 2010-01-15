### the directory name
set(directory source/VISUAL/VISUALIZER)

### list all filenames of the directory here
set(sources_list
AcquisitionInfoVisualizer.C
AcquisitionVisualizer.C
BaseVisualizer.C
BaseVisualizerGUI.C
ContactPersonVisualizer.C
DataProcessingVisualizer.C
DigestionVisualizer.C
DocumentIdentifierVisualizer.C
ExperimentalSettingsVisualizer.C
GradientVisualizer.C
HPLCVisualizer.C
InstrumentSettingsVisualizer.C
InstrumentVisualizer.C
IonDetectorVisualizer.C
IonSourceVisualizer.C
MassAnalyzerVisualizer.C
MetaInfoDescriptionVisualizer.C
MetaInfoVisualizer.C
ModificationVisualizer.C
PeptideHitVisualizer.C
PeptideIdentificationVisualizer.C
PrecursorVisualizer.C
ProteinHitVisualizer.C
ProteinIdentificationVisualizer.C
SampleVisualizer.C
ScanWindowVisualizer.C
SoftwareVisualizer.C
SourceFileVisualizer.C
SpectrumSettingsVisualizer.C
TaggingVisualizer.C
ProductVisualizer.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL\\VISUALIZER" FILES ${sources})

