### the directory name
set(directory source/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN)

### list all filenames of the directory here
set(sources_list
BackgroundControl.cpp
BackgroundIntensityBin.cpp
CentroidData.cpp
CentroidPeak.cpp
ClusteredMS2ConsensusSpectrum.cpp
Deisotoper.cpp
FTPeakDetectController.cpp
IsotopicDist.cpp
LCMSCData.cpp
LCMS.cpp
SuperHirnParameters.cpp
LCElutionPeak.cpp
MS1FeatureMerger.cpp
MS2ConsensusSpectrum.cpp
MS2Fragment.cpp
MS2Feature.cpp
ProcessData.cpp
RawData.cpp
ConsensusIsotopePattern.cpp
SHFeature.cpp
FeatureLCProfile.cpp
MS2Info.cpp
MSPeak.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(SuperHirn_sources ${SuperHirn_sources} ${sources} CACHE INTERNAL "This variable should hold all SuperHirn sources at the end of the config step" )

### source group definition
source_group("Source Files\\TRANSFORMATIONS\\FEATUREFINDER\\SUPERHIRN" FILES ${sources})
