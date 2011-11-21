### the directory name
set(directory source/FILTERING/DATAREDUCTION)

### list all filenames of the directory here
set(sources_list
DataFilters.C
ElutionPeakDetection.C
FeatureFindingMetabo.C
IsotopeDistributionCache.C
MassTraceDetection.C
SILACFilter.C
SILACFiltering.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FILTERING\\DATAREDUCTION" FILES ${sources})

