### the directory name
set(directory source/ANALYSIS/MAPMATCHING)

### list all filenames of the directory here
set(sources_list
BaseGroupFinder.C
BaseSuperimposer.C
ConsensusMapNormalizerAlgorithm.C
FeatureDistance.C
FeatureGroupingAlgorithm.C
FeatureGroupingAlgorithmIdentification.C
FeatureGroupingAlgorithmLabeled.C
FeatureGroupingAlgorithmQT.C
FeatureGroupingAlgorithmUnlabeled.C
LabeledPairFinder.C
MapAlignmentAlgorithm.C
MapAlignmentAlgorithmIdentification.C
MapAlignmentAlgorithmPoseClustering.C
MapAlignmentAlgorithmSpectrumAlignment.C
MapAlignmentEvaluationAlgorithm.C
MapAlignmentEvaluationAlgorithmPrecision.C
MapAlignmentEvaluationAlgorithmRecall.C
PoseClusteringAffineSuperimposer.C
PoseClusteringShiftSuperimposer.C
QTClusterFinder.C
SimplePairFinder.C
StablePairFinder.C
TransformationDescription.C
TransformationModel.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\MAPMATCHING" FILES ${sources})

