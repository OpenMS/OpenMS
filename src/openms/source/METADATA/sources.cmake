### the directory name
set(directory source/METADATA)

### list all filenames of the directory here
set(sources_list
Acquisition.cpp
AcquisitionInfo.cpp
ContactPerson.cpp
DataProcessing.cpp
Digestion.cpp
DocumentIdentifier.cpp
ExperimentalSettings.cpp
Gradient.cpp
HPLC.cpp
DocumentIDTagger.cpp
Instrument.cpp
InstrumentSettings.cpp
IonDetector.cpp
IonSource.cpp
MassAnalyzer.cpp
MetaInfo.cpp
MetaInfoDescription.cpp
MetaInfoInterface.cpp
MetaInfoRegistry.cpp
Modification.cpp
PeptideEvidence.cpp
PeptideHit.cpp
PeptideIdentification.cpp
Precursor.cpp
ProteinHit.cpp
ProteinIdentification.cpp
Sample.cpp
SampleTreatment.cpp
Software.cpp
SourceFile.cpp
SpectrumSettings.cpp
Tagging.cpp
ScanWindow.cpp
Product.cpp
Identification.cpp
SpectrumIdentification.cpp
IdentificationHit.cpp
ChromatogramSettings.cpp
CVTerm.cpp
CVTermList.cpp
MSQuantifications.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\METADATA" FILES ${sources})

