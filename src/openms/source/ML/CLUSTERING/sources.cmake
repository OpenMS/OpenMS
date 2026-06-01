### the directory name
set(directory source/ML/CLUSTERING)

### list all filenames of the directory here
set(sources_list
AverageLinkage.cpp
ClusterAnalyzer.cpp
ClusterFunctor.cpp
ClusterHierarchical.cpp
ClusteringGrid.cpp
CompleteLinkage.cpp
EuclideanSimilarity.cpp
GridBasedCluster.cpp
GridBasedClustering.cpp
SingleLinkage.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ML\\CLUSTERING" FILES ${sources})

