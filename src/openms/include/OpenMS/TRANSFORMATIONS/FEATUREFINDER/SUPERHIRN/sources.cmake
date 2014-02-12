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
FTPeakDetectController.h
IsotopicDist.h
LCMSCData.h
LCMS.h
SuperHirnParameters.h
LCElutionPeak.h
MS1FeatureMerger.h
MS2ConsensusSpectrum.h
MS2Fragment.h
MS2Feature.h
ProcessData.h
RawData.h
ConsensusIsotopePattern.h
SHFeature.h
FeatureLCProfile.h
MS2Info.h
MSPeak.h
SuperHirnUtil.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\TRANSFORMATIONS\\FEATUREFINDER\\SUPERHIRN" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

