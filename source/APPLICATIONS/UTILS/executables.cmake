### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(executables_list
CaapConvert
CaapEval
CVInspector
Digestor
FFEval
FuzzyDiff
HistView
IDExtractor
IdXMLInfo
LabeledEval
MapSimulator
RTEvaluation
SemanticValidator
SequenceCoverageCalculator
XMLValidator
IdXMLEvaluation
)

### pass source file list to the upper instance
set(UTILS_executables ${UTILS_executables} ${executables_list})

