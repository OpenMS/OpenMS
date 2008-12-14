### the directory name
set(directory source/APPLICATIONS/TOPP)

### list all filenames of the directory here
set(sources_list
AdditiveSeries.C
BaselineFilter.C
ConsensusID.C
DBExporter.C
DBImporter.C
DTAExtractor.C
Decharger.C
FalseDiscoveryRate.C
FeatureFinder.C
FeatureFinderMRM.C
FeatureLinker.C
FileConverter.C
FileFilter.C
FileInfo.C
FileMerger.C
IDDecoyProbability.C
IDFilter.C
IDMerger.C
INIFileEditor.C
ITRAQAnalyzer.C
InspectAdapter.C
InternalCalibration.C
MapAligner.C
MapNormalizer.C
MascotAdapter.C
MascotAdapterOnline.C
NoiseFilter.C
OMSSAAdapter.C
PILISIdentification.C
PILISModel.C
PTModel.C
PTPredict.C
PeakPicker.C
PepNovoAdapter.C
RTModel.C
RTPredict.C
Resampler.C
SILACAnalyzer.C
SequestAdapter.C
SpectraFilter.C
TOFCalibration.C
TOPPView.C
TextExporter.C
XTandemAdapter.C
splash.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### source group definition
source_group("Source Files\\APPLICATIONS\\TOPP" FILES ${sources})

