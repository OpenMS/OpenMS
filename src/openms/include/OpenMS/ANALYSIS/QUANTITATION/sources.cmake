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
ItraqConstants.h
ItraqEightPlexQuantitationMethod.h
ItraqFourPlexQuantitationMethod.h
KDTreeFeatureMaps.h
KDTreeFeatureNode.h
PeptideAndProteinQuant.h
ProteinInference.h
ProteinResolver.h
QuantitativeExperimentalDesign.h
TMTElevenPlexQuantitationMethod.h
TMTSixPlexQuantitationMethod.h
TMTTenPlexQuantitationMethod.h
AbsoluteQuantitation.h
AbsoluteQuantitationMethod.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\QUANTITATION" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
