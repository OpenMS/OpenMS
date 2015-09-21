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
FeatureGroupingAlgorithmIdentification.h
FeatureGroupingAlgorithmLabeled.h
FeatureGroupingAlgorithmQT.h
FeatureGroupingAlgorithmUnlabeled.h
LabeledPairFinder.h
MapAlignmentAlgorithmIdentification.h
MapAlignmentAlgorithmPoseClustering.h
MapAlignmentAlgorithmSpectrumAlignment.h
MapAlignmentEvaluationAlgorithm.h
MapAlignmentEvaluationAlgorithmPrecision.h
MapAlignmentEvaluationAlgorithmRecall.h
MapAlignmentTransformer.h
PoseClusteringAffineSuperimposer.h
PoseClusteringShiftSuperimposer.h
QTClusterFinder.h
SimplePairFinder.h
StablePairFinder.h
TransformationDescription.h
TransformationModel.h
TransformationModelBSpline.h
TransformationModelLinear.h
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

