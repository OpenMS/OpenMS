### the directory name
set(directory source/ANALYSIS/QUANTITATION)

### list all filenames of the directory here
set(sources_list
IsobaricChannelExtractor.cpp
IsobaricIsotopeCorrector.cpp
IsobaricNormalizer.cpp
IsobaricQuantifier.cpp
IsobaricQuantifierStatistics.cpp
IsobaricQuantitationMethod.cpp
ItraqChannelExtractor.cpp
ItraqConstants.cpp
ItraqEightPlexQuantitationMethod.cpp
ItraqFourPlexQuantitationMethod.cpp
ItraqQuantifier.cpp
PeptideAndProteinQuant.cpp
ProteinInference.cpp
ProteinResolver.cpp
QuantitativeExperimentalDesign.cpp
TMTSixPlexQuantitationMethod.cpp
TMTTenPlexQuantitationMethod.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\QUANTITATION" FILES ${sources})
