### the directory name
set(directory source/ANALYSIS/TARGETED)

### list all filenames of the directory here
set(sources_list
MRMMapping.cpp
PrecursorIonSelectionPreprocessing.cpp
PrecursorIonSelection.cpp
OfflinePrecursorIonSelection.cpp
PSLPFormulation.cpp
PSProteinInference.cpp
IncludeExcludeTarget.cpp
TargetedExperiment.cpp
TargetedExperimentHelper.cpp
InclusionExclusionList.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\TARGETED" FILES ${sources})

