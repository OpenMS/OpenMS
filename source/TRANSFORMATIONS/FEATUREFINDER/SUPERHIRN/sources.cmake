### the directory name
set(directory source/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN)

### list all filenames of the directory here
set(sources_list
BackgroundControl.C
BackgroundIntensityBin.C
CentroidData.C
CentroidPeak.C
ClusteredMS2ConsensusSpectrum.C
Deisotoper.C
ExternalIsotopicDistribution.C
FT_PEAK_DETEC_mzXML_reader.C
FT_PeakDetectController.C
IsotopicDist.C
LCMSCData.C
LC_MS.C
LC_MS_XML_reader.C
LC_elution_peak.C
MS1_feature_merger.C
MS2ConsensusSpectrum.C
MS2Fragment.C
MS2_feature.C
PeptideIsotopeDistribution.C
Process_Data.C
RawData.C
consensIsotopePattern.C
feature.C
featureLCprofile.C
ms2_info.C
ms_peak.C
simple_math2.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\TRANSFORMATIONS\\FEATUREFINDER\\SUPERHIRN" FILES ${sources})
