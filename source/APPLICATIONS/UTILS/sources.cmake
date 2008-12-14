### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(sources_list
CAAP_groundtruth_to_consensusXML.C
CVInspector.C
Digestor.C
FFEval.C
FuzzyDiff.C
HistView.C
IDExtractor.C
IdXMLInfo.C
LabeledEval.C
RTEvaluation.C
SemanticValidator.C
SequenceCoverageCalculator.C
XMLValidator.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\APPLICATIONS\\UTILS" FILES ${sources})

