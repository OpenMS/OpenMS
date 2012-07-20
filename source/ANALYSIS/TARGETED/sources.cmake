### the directory name
set(directory source/ANALYSIS/TARGETED)

### list all filenames of the directory here
set(sources_list
PrecursorIonSelectionPreprocessing.C
PrecursorIonSelection.C
OfflinePrecursorIonSelection.C
PSLPFormulation.C
PSProteinInference.C
IncludeExcludeTarget.C
TargetedExperiment.C
TargetedExperimentHelper.C
InclusionExclusionList.C
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

