### the directory name
set(directory include/OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN)

### list all header files of the directory here
set(sources_list_h
BackgroundControl.h
BackgroundIntensityBin.h
CentroidData.h
CentroidPeak.h
ClusteredMS2ConsensusSpectrum.h
Deisotoper.h
ExternalIsotopicDistribution.h
FT_PEAK_DETEC_mzXML_reader.h
FT_PeakDetectController.h
IsotopicDist.h
LCMSCData.h
LC_MS.h
LC_MS_XML_reader.h
LC_elution_peak.h
MS1_feature_merger.h
MS2ConsensusSpectrum.h
MS2Fragment.h
MS2_feature.h
PeptideIsotopeDistribution.h
Process_Data.h
RawData.h
consensIsotopePattern.h
feature.h
featureLCprofile.h
ms2_info.h
ms_peak.h
simple_math2.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\TRANSFORMATIONS\\FEATUREFINDER\\SUPERHIRN" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

