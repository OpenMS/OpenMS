### the directory name
set(directory include/OpenMS/ANALYSIS/MAPMATCHING)

### list all header files of the directory here
set(sources_list_h
BaseGroupFinder.h
BaseSuperimposer.h
ConsensusMapNormalizerAlgorithmThreshold.h
ConsensusMapNormalizerAlgorithmMedian.h
ConsensusMapNormalizerAlgorithmQuantile.h
FeatureDistance.h
FeatureGroupingAlgorithm.h
FeatureGroupingAlgorithmLabeled.h
FeatureGroupingAlgorithmKD.h
FeatureGroupingAlgorithmQT.h
FeatureGroupingAlgorithmUnlabeled.h
FeatureMapping.h
LabeledPairFinder.h
MapAlignmentAlgorithmIdentification.h
MapAlignmentAlgorithmKD.h
MapAlignmentAlgorithmPoseClustering.h
MapAlignmentAlgorithmTreeGuided.h
MapAlignmentEvaluationAlgorithm.h
MapAlignmentEvaluationAlgorithmPrecision.h
MapAlignmentEvaluationAlgorithmRecall.h
MapAlignmentTransformer.h
PoseClusteringAffineSuperimposer.h
PoseClusteringShiftSuperimposer.h
QTClusterFinder.h
StablePairFinder.h
TransformationDescription.h
TransformationModel.h
TransformationModelBSpline.h
TransformationModelLinear.h
TransformationModelLowess.h
TransformationModelInterpolated.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\MAPMATCHING" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

