### the directory name
set(directory include/OpenMS/DATASTRUCTURES)

### list all header files of the directory here
set(sources_list_h
Adduct.h
BinaryTreeNode.h
CalibrationData.h
ChargePair.h
Compomer.h
ConstRefVector.h
ConvexHull2D.h
CVMappingTerm.h
CVMappingRule.h
CVReference.h
CVMappings.h
DBoundingBox.h
DIntervalBase.h
DPosition.h
DRange.h
DataValue.h
Date.h
DateTime.h
DefaultParamHandler.h
DistanceMatrix.h
FASTAContainer.h
GridFeature.h
IsotopeCluster.h
KDTree.h
ListUtils.h
ListUtilsIO.h
LPWrapper.h
Map.h
MassExplainer.h
MatchedIterator.h
Matrix.h
Param.h
QTCluster.h
SeqanIncludeWrapper.h
String.h
StringUtils.h
StringListUtils.h
ToolDescription.h

)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\DATASTRUCTURES" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
