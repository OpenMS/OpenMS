### the directory name
set(directory include/OpenMS/DATASTRUCTURES)

### list all header files of the directory here
set(sources_list_h
Adduct.h
BigString.h
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
DoubleList.h
GridElement.h
GridFeature.h
HashGrid.h
IntList.h
IsotopeCluster.h
Map.h
MassExplainer.h
Matrix.h
Param.h
SeqanIncludeWrapper.h
SparseVector.h
String.h
StringList.h
SuffixArray.h
SuffixArrayPeptideFinder.h
SuffixArraySeqan.h
SuffixArrayTrypticCompressed.h
SuffixArrayTrypticSeqan.h
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

