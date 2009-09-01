### the directory name
set(directory include/OpenMS/METADATA)

### list all header files of the directory here
set(sources_list_h
Acquisition.h
AcquisitionInfo.h
ContactPerson.h
DataProcessing.h
Digestion.h
DocumentIdentifier.h
ExperimentalSettings.h
Gradient.h
HPLC.h
IDTagger.h
Instrument.h
InstrumentSettings.h
IonDetector.h
IonSource.h
MassAnalyzer.h
MetaInfo.h
MetaInfoDescription.h
MetaInfoInterface.h
MetaInfoRegistry.h
Modification.h
PeptideHit.h
PeptideIdentification.h
Precursor.h
ProteinHit.h
ProteinIdentification.h
Sample.h
SampleTreatment.h
Software.h
SourceFile.h
SpectrumSettings.h
Tagging.h
ScanWindow.h
Product.h
PeptideEvidence.h
Identification.h
SpectrumIdentification.h
IdentificationHit.h
ChromatogramSettings.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\METADATA" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

