### the directory name
set(directory source/ANALYSIS/MAPMATCHING)

### list all filenames of the directory here
set(sources_list
BaseGroupFinder.cpp
BaseSuperimposer.cpp
ConsensusMapNormalizerAlgorithmThreshold.cpp
ConsensusMapNormalizerAlgorithmMedian.cpp
ConsensusMapNormalizerAlgorithmQuantile.cpp
FeatureDistance.cpp
FeatureGroupingAlgorithm.cpp
FeatureGroupingAlgorithmIdentification.cpp
FeatureGroupingAlgorithmLabeled.cpp
FeatureGroupingAlgorithmQT.cpp
FeatureGroupingAlgorithmUnlabeled.cpp
LabeledPairFinder.cpp
MapAlignmentAlgorithmIdentification.cpp
MapAlignmentAlgorithmPoseClustering.cpp
MapAlignmentAlgorithmSpectrumAlignment.cpp
MapAlignmentEvaluationAlgorithm.cpp
MapAlignmentEvaluationAlgorithmPrecision.cpp
MapAlignmentEvaluationAlgorithmRecall.cpp
MapAlignmentTransformer.cpp
PoseClusteringAffineSuperimposer.cpp
PoseClusteringShiftSuperimposer.cpp
QTClusterFinder.cpp
SimplePairFinder.cpp
StablePairFinder.cpp
TransformationDescription.cpp
TransformationModel.cpp
TransformationModelBSpline.cpp
TransformationModelLinear.cpp
TransformationModelInterpolated.cpp
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

