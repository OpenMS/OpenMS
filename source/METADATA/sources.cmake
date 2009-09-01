### the directory name
set(directory source/METADATA)

### list all filenames of the directory here
set(sources_list
Acquisition.C
AcquisitionInfo.C
ContactPerson.C
DataProcessing.C
Digestion.C
DocumentIdentifier.C
ExperimentalSettings.C
Gradient.C
HPLC.C
IDTagger.C
Instrument.C
InstrumentSettings.C
IonDetector.C
IonSource.C
MassAnalyzer.C
MetaInfo.C
MetaInfoDescription.C
MetaInfoInterface.C
MetaInfoRegistry.C
Modification.C
PeptideHit.C
PeptideIdentification.C
Precursor.C
ProteinHit.C
ProteinIdentification.C
Sample.C
SampleTreatment.C
Software.C
SourceFile.C
SpectrumSettings.C
Tagging.C
ScanWindow.C
Product.C
PeptideEvidence.C
Identification.C
SpectrumIdentification.C
IdentificationHit.C
ChromatogramSettings.C
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

