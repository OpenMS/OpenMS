### the directory name
set(directory include/OpenMS/COMPARISON/CLUSTERING)

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
HashGrid.h
SingleLinkage.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\COMPARISON\\CLUSTERING" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

