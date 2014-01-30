### the directory name
set(directory source/ANALYSIS/QUANTITATION)

### list all filenames of the directory here
set(sources_list
ItraqChannelExtractor.cpp
ItraqConstants.cpp
ItraqQuantifier.cpp
PeptideAndProteinQuant.cpp
ProteinInference.cpp
ProteinResolver.cpp
QuantitativeExperimentalDesign.cpp
IsobaricQuantitationMethod.cpp
IsobaricChannelExtractor.cpp
ItraqFourPlexQuantitationMethod.cpp
IsobaricQuantifier.cpp
IsobaricNormalizer.cpp
IsobaricQuantifierStatistics.cpp
IsobaricIsotopeCorrector.cpp
ItraqEightPlexQuantitationMethod.cpp
TMTSixPlexQuantitationMethod.cpp
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

