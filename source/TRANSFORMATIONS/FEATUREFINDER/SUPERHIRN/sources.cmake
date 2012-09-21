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
FTPeakDetectController.C
IsotopicDist.C
LCMSCData.C
LCMS.C
SuperHirnParameters.C
LCElutionPeak.C
MS1FeatureMerger.C
MS2ConsensusSpectrum.C
MS2Fragment.C
MS2Feature.C
ProcessData.C
RawData.C
ConsensusIsotopePattern.C
SHFeature.C
FeatureLCProfile.C
MS2Info.C
MSPeak.C
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
