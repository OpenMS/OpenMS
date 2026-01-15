### the directory name
set(directory include/OpenMS/MATH/STATISTICS)

### list all header files of the directory here
set(sources_list_h
BasicStatistics.h
GammaDistributionFitter.h
GaussFitter.h
GumbelDistributionFitter.h
GumbelMaxLikelihoodFitter.h
Histogram.h
PosteriorErrorProbabilityModel.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\MATH\\STATISTICS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
