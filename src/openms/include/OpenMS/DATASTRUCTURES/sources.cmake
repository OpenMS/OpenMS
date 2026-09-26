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
ExposedVector.h
FASTAContainer.h
FlagSet.h
GridFeature.h
IsotopeCluster.h
KDTree.h
ListUtils.h
ListUtilsIO.h
LPWrapper.h
MassExplainer.h
MatchedIterator.h
Matrix.h
OSWData.h
Param.h
ParamTags.h
ParamValue.h
QTCluster.h
StringUtils.h
StringListUtils.h
ToolDescription.h
ToolInfo.h
TypeAliases.h
RegularExpression.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\DATASTRUCTURES" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

### Private (non-installed) header: MatrixEigen.h includes <Eigen/Core> and its
### eigenView()/eigenVectorView()/eigenMatrixView() helpers return Eigen::Map
### types. It says so itself ("INTERNAL header file ... should NOT be included in
### public headers"), and only libOpenMS sources and class tests include it.
### Keeping it off OpenMS_sources_h is what lets Eigen3::Eigen be a PRIVATE link
### dependency: no Eigen type or include appears in any installed header.
set(private_headers_list_h
MatrixEigen.h
)
set(private_sources_h)
foreach(i ${private_headers_list_h})
	list(APPEND private_sources_h ${directory}/${i})
endforeach(i)
source_group("Header Files\\OpenMS\\DATASTRUCTURES" FILES ${private_sources_h})
set(OpenMS_private_headers ${OpenMS_private_headers} ${private_sources_h})
