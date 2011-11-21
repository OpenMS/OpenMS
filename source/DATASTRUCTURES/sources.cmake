### the directory name
set(directory source/DATASTRUCTURES)

### list all filenames of the directory here
set(sources_list
Adduct.C
BigString.C
BinaryTreeNode.C
ChargePair.C
Compomer.C
ConstRefVector.C
ConvexHull2D.C
CVMappingTerm.C
CVMappingRule.C
CVReference.C
CVMappings.C
DBoundingBox.C
DIntervalBase.C
DPosition.C
DRange.C
DataValue.C
Date.C
DateTime.C
DefaultParamHandler.C
DistanceMatrix.C
DoubleList.C
GridFeature.C
IntList.C
Map.C
MassExplainer.C
Matrix.C
Param.C
QTCluster.C
SparseVector.C
String.C
StringList.C
SuffixArray.C
SuffixArrayPeptideFinder.C
SuffixArraySeqan.C
SuffixArrayTrypticCompressed.C
SuffixArrayTrypticSeqan.C
ToolDescription.C
LPWrapper.C
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

