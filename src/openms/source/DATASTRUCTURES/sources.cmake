### the directory name
set(directory source/DATASTRUCTURES)

### list all filenames of the directory here
set(sources_list
Adduct.cpp
BigString.cpp
BinaryTreeNode.cpp
CalibrationData.cpp
ChargePair.cpp
Compomer.cpp
ConstRefVector.cpp
ConvexHull2D.cpp
CVMappingTerm.cpp
CVMappingRule.cpp
CVReference.cpp
CVMappings.cpp
DBoundingBox.cpp
DIntervalBase.cpp
DPosition.cpp
DRange.cpp
DataValue.cpp
Date.cpp
DateTime.cpp
DefaultParamHandler.cpp
DistanceMatrix.cpp
FASTAContainer.cpp
GridFeature.cpp
ListUtils.cpp
ListUtilsIO.cpp
LPWrapper.cpp
Map.cpp
MassExplainer.cpp
Matrix.cpp
Param.cpp
QTCluster.cpp
SparseVector.cpp
String.cpp
StringListUtils.cpp
StringUtils.cpp
ToolDescription.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\DATASTRUCTURES" FILES ${sources})

