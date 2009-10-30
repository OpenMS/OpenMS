### the directory name
set(directory source/APPLICATIONS/TOPP)

### list all filenames of the directory here
set(TOPP_executables
AdditiveSeries
BaselineFilter
ConsensusID
DBExporter
DBImporter
DTAExtractor
Decharger
FalseDiscoveryRate
FeatureFinder
FeatureLinker
FileConverter
FileFilter
FileInfo
FileMerger
IDDecoyProbability
IDFilter
IDMerger
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
PTModel
PTPredict
PeakPicker
PepNovoAdapter
RTModel
RTPredict
Resampler
PILISModel
PILISIdentification
SequestAdapter
SpecLibSearcher
SpectraFilter
TOFCalibration
TOPPView
TOPPAS
TextExporter
TextImporter
XTandemAdapter
IDMapper
IDRTCalibration
SILACAnalyzer
PrecursorIonSelector
IDFileConverter
CompNovo
PeptideIndexer
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TOPP_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})
