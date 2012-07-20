### the directory name
set(directory include/OpenMS/ANALYSIS/TARGETED)

### list all header files of the directory here
set(sources_list_h
PrecursorIonSelectionPreprocessing.h
PrecursorIonSelection.h
OfflinePrecursorIonSelection.h
PSLPFormulation.h
PSProteinInference.h
IncludeExcludeTarget.h
TargetedExperiment.h
TargetedExperimentHelper.h
InclusionExclusionList.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\TARGETED" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

