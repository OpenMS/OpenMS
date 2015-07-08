### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(UTILS_executables
AccurateMassSearch
CVInspector
DecoyDatabase
DeMeanderize
Digestor
DigestorMotif
ERPairFinder
FeatureFinderSuperHirn
FFEval
FuzzyDiff
IDDecoyProbability
IDExtractor
IDMassAccuracy
IDScoreSwitcher
IDSplitter
LabeledEval
LowMemPeakPickerHiRes
LowMemPeakPickerHiRes_RandomAccess
MapAlignmentEvaluation
MassCalculator
MetaboliteSpectralMatcher
MetaProSIP
MRMPairFinder
MSSimulator
OpenMSInfo
PeakPickerIterative
QCCalculator
QCEmbedder
QCExporter
QCExtractor
QCImporter
QCMerger
QCShrinker
RNPxl
RNPxlXICFilter
RTAnnotator
RTEvaluation
SemanticValidator
SequenceCoverageCalculator
SimpleSearchEngine
SpecLibCreator
SvmTheoreticalSpectrumGeneratorTrainer
TransformationEvaluation
TopPerc
XMLValidator
#SimpleSearchEngine
#RNPxlSearch
)

if(NOT DISABLE_OPENSWATH)
  set(UTILS_executables
    ${UTILS_executables}
    ConvertTSVToTraML
    ConvertTraMLToTSV
    OpenSwathDIAPreScoring
    OpenSwathMzMLFileCacher
    OpenSwathWorkflow
    OpenSwathRewriteToFeatureXML
    MRMTransitionGroupPicker
  )
endif(NOT DISABLE_OPENSWATH)


## all targets requiring OpenMS_GUI
set(UTILS_executables_with_GUIlib
IDEvaluator
ImageCreator
INIUpdater
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${UTILS_executables} ${UTILS_executables_with_GUIlib})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)

source_group("" FILES ${sources_VS})
