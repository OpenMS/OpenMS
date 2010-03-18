### the directory name
set(directory source/APPLICATIONS/TOPP)

### list all filenames of the directory here
set(TOPP_executables
AdditiveSeries
BaselineFilter
CompNovo
ConsensusID
DBExporter
DBImporter
DTAExtractor
Decharger
ExecutePipeline
FalseDiscoveryRate
FeatureFinder
FeatureLinker
FileConverter
FileFilter
FileInfo
FileMerger
GenericWrapper
IDDecoyProbability
IDPosteriorErrorProbability
IDFileConverter
IDFilter
IDMapper
IDMerger
IDRTCalibration
INIFileEditor
ITRAQAnalyzer
InspectAdapter
InternalCalibration
MapAligner
MapNormalizer
MascotAdapter
MascotAdapterOnline
NoiseFilter
OMSSAAdapter
PILISIdentification
PILISModel
PTModel
PTPredict
PeakPicker
PepNovoAdapter
PeptideIndexer
PrecursorIonSelector
PrecursorMassCorrector
ProteinInference
ProteinQuantifier
RTModel
RTPredict
Resampler
SILACAnalyzer
SILACAnalyzer2
SeedListGenerator
SequestAdapter
SpecLibSearcher
SpectraFilter
TOFCalibration
TOPPAS
TOPPView
TextExporter
XTandemAdapter
InclusionExclusionListCreator
SpectraMerger
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TOPP_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})
