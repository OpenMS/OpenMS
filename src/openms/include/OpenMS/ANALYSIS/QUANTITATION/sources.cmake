### the directory name
set(directory include/OpenMS/ANALYSIS/QUANTITATION)

### list all header files of the directory here
set(sources_list_h
IsobaricChannelExtractor.h
IsobaricIsotopeCorrector.h
IsobaricNormalizer.h
IsobaricQuantifier.h
IsobaricQuantifierStatistics.h
IsobaricQuantitationMethod.h
ItraqChannelExtractor.h
ItraqConstants.h
ItraqEightPlexQuantitationMethod.h
ItraqFourPlexQuantitationMethod.h
ItraqQuantifier.h
PeptideAndProteinQuant.h
ProteinInference.h
ProteinResolver.h
QuantitativeExperimentalDesign.h
TMTSixPlexQuantitationMethod.h
TMTTenPlexQuantitationMethod.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\QUANTITATION" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
