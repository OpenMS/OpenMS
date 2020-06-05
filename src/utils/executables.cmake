### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(UTILS_executables
AccurateMassSearch
AssayGeneratorMetabo
ClusterMassTraces
ClusterMassTracesByPrecursor
CVInspector
DatabaseFilter
DecoyDatabase
DeMeanderize
Digestor
DigestorMotif
Epifany
ERPairFinder
FeatureFinderMetaboIdent
FeatureFinderSuperHirn
FFEval
FuzzyDiff
IDDecoyProbability
IDExtractor
IDMassAccuracy
IDScoreSwitcher
IDSplitter
LabeledEval
MassCalculator
MetaboliteAdductDecharger
MetaboliteSpectralMatcher
MetaProSIP
MRMPairFinder
MSFraggerAdapter
MSSimulator
MSstatsConverter
MultiplexResolver
MzMLSplitter
NovorAdapter
NucleicAcidSearchEngine
OpenMSInfo
PeakPickerIterative
PSMFeatureExtractor
QCCalculator
QCEmbedder
QCExporter
QCExtractor
QCImporter
QCMerger
QCShrinker
ProteomicsLFQ
RNADigestor
RNAMassCalculator
RNPxlXICFilter
RNPxlSearch
RTEvaluation
SemanticValidator
SequenceCoverageCalculator
SimpleSearchEngine
SiriusAdapter
SpecLibCreator
SpectraSTSearchAdapter
SvmTheoreticalSpectrumGeneratorTrainer
TICCalculator
TransformationEvaluation
XMLValidator
)

if(NOT DISABLE_OPENSWATH)
  set(UTILS_executables
    ${UTILS_executables}
    TargetedFileConverter
    OpenSwathDIAPreScoring
    OpenSwathMzMLFileCacher
    OpenSwathWorkflow
    OpenSwathFileSplitter
    OpenSwathRewriteToFeatureXML
    MRMTransitionGroupPicker
  )
endif(NOT DISABLE_OPENSWATH)


## all targets requiring OpenMS_GUI
set(UTILS_executables_with_GUIlib
ImageCreator
INIUpdater
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${UTILS_executables} ${UTILS_executables_with_GUIlib})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)

source_group("" FILES ${sources_VS})
