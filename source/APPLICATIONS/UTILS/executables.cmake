### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(UTILS_executables
ConvertTSVToTraML
ConvertTraMLToTSV
CVInspector
DeMeanderize
DecoyDatabase
Digestor
DigestorMotif
ERPairFinder
FeatureFinderSuperHirn
FFEval
FuzzyDiff
IDEvaluation
IDExtractor
IDMassAccuracy
IDSplitter
IDDecoyProbability
RTEvaluation
ImageCreator
INIUpdater
LabeledEval
MassCalculator
MRMPairFinder
MRMMapper
MSSimulator
MapAlignmentEvaluation
OpenMSInfo
OpenSwathDecoyGenerator
OpenSwathFeatureXMLToTSV
OpenSwathChromatogramExtractor
OpenSwathAnalyzer
OpenSwathRTNormalizer
OpenSwathmzMLFileCacher
SemanticValidator
SequenceCoverageCalculator
SpecLibCreator
SvmTheoreticalSpectrumGeneratorTrainer
TransformationEvaluation
XMLValidator
QCCalculator
)

## all targets with need linkage against OpenMS_GUI.lib - they also need to appear in the list above)
set(UTILS_executables_with_GUIlib
IDEvaluation
ImageCreator
INIUpdater
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${UTILS_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})
