### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(UTILS_executables
AccurateMassSearch
CVInspector
DeMeanderize
DecoyDatabase
Digestor
DigestorMotif
ERPairFinder
FeatureFinderSuperHirn
FFEval
FuzzyDiff
IDEvaluator
IDExtractor
IDMassAccuracy
IDSplitter
IDDecoyProbability
RTEvaluation
ImageCreator
INIUpdater
LabeledEval
LowMemPeakPickerHiRes
MassCalculator
MRMPairFinder
MSSimulator
MapAlignmentEvaluation
OpenMSInfo
PeakPickerRapid
PeakPickerIterative
SemanticValidator
SequenceCoverageCalculator
SpecLibCreator
SvmTheoreticalSpectrumGeneratorTrainer
TransformationEvaluation
XMLValidator
QCCalculator
QCImporter
QCEmbedder
QCExporter
QCExtractor
QCMerger
QCShrinker
RNPxl
RNPxlXICFilter
)

if(NOT DISABLE_OPENSWATH)
  set(UTILS_executables
    ${UTILS_executables}
    ConvertTSVToTraML
    ConvertTraMLToTSV
    OpenSwathDIAPreScoring
    OpenSwathMzMLFileCacher
    OpenSwathRewriteToFeatureXML
    MRMTransitionGroupPicker
  )
endif(NOT DISABLE_OPENSWATH)


## all targets with need linkage against OpenMS_GUI.lib - they also need to appear in the list above)
set(UTILS_executables_with_GUIlib
IDEvaluator
ImageCreator
INIUpdater
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${UTILS_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})
