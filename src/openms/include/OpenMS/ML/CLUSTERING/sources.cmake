### the directory name
set(directory include/OpenMS/ML/CLUSTERING)

### list all header files of the directory here
set(sources_list_h
AverageLinkage.h
ClusterAnalyzer.h
ClusterFunctor.h
ClusterHierarchical.h
ClusteringGrid.h
CompleteLinkage.h
EuclideanSimilarity.h
GridBasedCluster.h
GridBasedClustering.h
SingleLinkage.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ML\\CLUSTERING" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

### Private (non-installed) header: HashGrid uses Boost.Unordered because QTClusterFinder
### derives its cluster order from the grid's bucket order, and only Boost gives the same
### order on every platform (the std containers' bucket policies differ, which showed up as
### a Windows-only FeatureLinkerUnlabeledQT regression). QTClusterFinder hands out an opaque
### Grid instead of this template, so Boost stays out of the installed interface.
set(private_sources_h ${directory}/HashGrid.h)
source_group("Header Files\\OpenMS\\ML\\CLUSTERING" FILES ${private_sources_h})
set(OpenMS_private_headers ${OpenMS_private_headers} ${private_sources_h})

