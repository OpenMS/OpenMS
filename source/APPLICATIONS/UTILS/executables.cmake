### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(executables_list
CaapConvert
CVInspector
Digestor
FFEval
FuzzyDiff
HistView
IDExtractor
LabeledEval
MapAlignmentEvaluation
MSSimulator
RTEvaluation
SemanticValidator
SequenceCoverageCalculator
XMLValidator
IdXMLEvaluation
DecoyDatabase
IDMassAccuracy
PeptideIndexer
)

### pass source file list to the upper instance
set(UTILS_executables ${UTILS_executables} ${executables_list})

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${executables_list})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})
